#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 05-12-2025
#
#   Optimizes DADA2 in two ways:
#       1. Splits mapping file in batches, can be specified via --batch_n
#       2. Disables multithreading at denoise but allocates each sample on a separate CPU
#         2a. Option to change to classic way is still possible by -p, --parallel flag
#
#   Required INPUT: mapping file
#
#   OUTPUT:
#     1. denoising_stats as RDS file
#     2. seq_tab as RDS file
#     3. ErrProfile
#
#   Requires Qiime2 to import the data in the right format
#
#   Containers:
#     - Dockerfile to build docker for parallel_dada2.R
#
#------------------------------------------------------------------------------------#

Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

# required libraries & Loess functions
library("foreach")
library("dplyr")

# Try both paths, docker or git 
paths_to_try <- c(paste0(current_path, "R/error_methods.R"), "error_methods.R")
for (path in paths_to_try) {
  if (file.exists(path)) {
    source(path)
    break
  }
}

#-----------------------------------------#
# Parsing from command line               #
#-----------------------------------------#
option_list <- list (optparse::make_option(c("-m", "--metadata"),
                                           action = "store",
                                           help="tab seperated file"),
                     optparse::make_option(c("-n", "--batch_n"),
                                           action = "store",
                                           default = 500,
                                           help="Specify batch size"),
                     optparse::make_option(c("-c", "--cpus"),
                                           action = "store",
                                           default = 8,
                                           help="Specify number of cores to be used"),
                     optparse::make_option(c("-s", "--seed"),
                                           action = "store",
                                           default = 100,
                                           help = "sets seed number"),
                     optparse::make_option(c("-p", "--parallel"),
                                           action = "store_true",
                                           default = FALSE,
                                           help = "parallel TRUE will allocate sample per core, FALSE will split sample among cores"),
                     optparse::make_option(c("--novaseq"),
                                           action = "store_true",
                                           default = FALSE,
                                           help = "Specify if sequence data originates from novaseq"),

                     # Optional arguments for dada2 from command line
                     optparse::make_option(c("--p-trunc-len"),
                                           action = "store",
                                           default = 0,
                                           help = "Default 0. Truncate reads after truncLen bases. Reads shorter than this are discarded."),
                     optparse::make_option(c("--p-trunc-q"),
                                           action = "store",
                                           default = 2,
                                           help = "Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ"),
                     optparse::make_option(c("--p-max-ee"),
                                           action = "store",
                                           default = Inf,
                                           help = "Default Inf (no EE filtering). After truncation, reads with higher than maxEE 'expected errors' will be discarded."),
                     optparse::make_option(c("--p-min-fold-parent-over-abundance"),
                                           action = "store",
                                           default = 1,
                                           help = "Values should be greater than or equal to 1 (i.e. parents should be more abundant than the sequence being tested)."),
                     optparse::make_option(c("--p-chimera-method"),
                                           action = "store",
                                           default = "consensus",
                                           help = "Default is 'consensus'. Only has an effect if a sequence table is provided. Options: 'pooled', 'consensus', 'per-sample' see dada2 docs"),
                     optparse::make_option(c("--skip-denoise"),
                                           action = "store_true",
                                           default = FALSE,
                                           help = "You would never skip denoise, unless the input data is not compatible to dada2. For example you have nanopore long-reads (Default: FALSE)")
)

# Collects arguments
parser <- optparse::OptionParser(option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments=TRUE)

opt <- arguments$options

#-----------------------------------------#
# Loads required files                    #
#-----------------------------------------#

if (!is.null(opt$metadata)) {
  mapping <- data.table::fread(opt$metadata, header = TRUE)
} else stop("Please provide a tab-separated metadata file!")

# Fetch user-input or default parameters
batch_size <- opt$batch_n
seed_n <- opt$seed
cpus_n <- opt$cpus

#-----------------------------------------#
# Required functions                      #
#-----------------------------------------#
# Counts N of reads
getN <- function(x) sum(dada2::getUniques(x))

add_percentage_columns <- function(df) {
  denom_col <- df[[1]]
  for (col in names(df)) {
    if (is.numeric(df[[col]])) {
      pct_col_name <- paste0(col, " [%]")
      df[[pct_col_name]] <- (df[[col]] / denom_col) * 100
    }
  }
  return(df)
}

#-----------------------------------------#
# Setting up parallel and seed            #
#-----------------------------------------#

# Setting up seed
set.seed(seed = seed_n)

if (opt$parallel) {
  # Setting up parallel backend
  cl <- parallel::makeCluster(cpus_n)
  doParallel::registerDoParallel(cl)
}

#-----------------------------------------#
# Preparation of batches                  #
#-----------------------------------------#

# Shuffling rows
mapping[base::sample(.N)]
colnames(mapping) <- tolower(colnames(mapping))
mapping_n <- nrow(mapping)

# Automatically returns a single batch if mapping_n < batch_size
batches <- base::split(mapping, base::ceiling(seq_len(mapping_n) / batch_size))

# Required paths
filtpath <- file.path("filtered")

#-----------------------------------------#
# Runs denoise in batches                 #
#-----------------------------------------#
for (i in 1:length(batches)) {
  # Assigning sample names and fastq path from mapping
  sample_names <- batches[[i]][["sample-id"]]
  sample_fastq <- batches[[i]][["absolute-filepath"]]

  # Setting filtered paths
  filtFs <- file.path(filtpath, basename(sample_fastq))

  # Filtering script
  out <- dada2::filterAndTrim(
    fwd = sample_fastq,
    filt = filtFs,
    truncLen = opt$`p-trunc-len`,
    maxEE = opt$`p-max-ee`,
    truncQ = opt$`p-trunc-q`,
    rm.phix = TRUE,
    compress = TRUE,
    verbose = TRUE,
    multithread = cpus_n
    )

  # Dereplication
  derep <- dada2::derepFastq(filtFs, verbose = TRUE)
  if (length(sample_names) == 1) {
    derepFs <- list(derep)
  } else {
    derepFs <- derep
  }
  rm(derep)
  names(derepFs) <- sample_names

  if (!opt$`skip-denoise`) {
    # Learn error rates
    if (opt$novaseq) {
      # I choose model 4 based on pre-liminary pilot tests of deeply sequenced NovaSeq data (800k - 3 million reads)
      err <- dada2::learnErrors(
        fls = derepFs,
        multithread = cpus_n,
        errorEstimationFunction = loessErrfun_mod4
        )
    } else {
      err <- dada2::learnErrors(
        fls = derepFs,
        multithread = cpus_n
        )
    }


    # Save err plot
    ggplot2::ggsave(filename = paste0("errProfile_", i, ".png"),
                    plot = dada2::plotErrors(err, nominalQ=TRUE),
                    width = 10,
                    height = 10,
                    dpi = 400)

    # Parallel Denoising
    if (opt$parallel) {
      dds <- foreach::foreach(sam = sample_names, .combine = "c", .packages = "dada2") %dopar% {
        cat("Processing:", sam, "\n")
        list(sam = dada2::dada(
          derep = derepFs[[ sam ]],
          err = err,
          multithread = FALSE
          )
        )
      }
    } else {
      dds <- vector("list", length(sample_names))
      names(dds) <- sample_names
      for (sam in sample_names) {
        dds[[sam]] <- dada2::dada(
          derep = derepFs[[ sam ]],
          err = err,
          multithread = cpus_n
          )
      }
    }

    # Create sequence table
    dds <- dds[!sapply(dds, is.null)]
    seqtab <- dada2::makeSequenceTable(dds)
    rownames(seqtab) <- sample_names

    # stats of reads
    track <- cbind(out, sapply(dds, getN))
    colnames(track) <- c("input", "filtered", "denoised")
    rownames(track) <- sample_names

  } else {

    ## Creating seqtab from dereplicated files
    features <- unique(unlist(lapply(derepFs, function(x) names(x$uniques))))

    ## Build otu table
    otu_tab <- sapply(derepFs, function(x) {
      counts <- x$uniques
      as.integer(counts[features])
    })
    otu_tab[is.na(otu_tab)] <- 0
    otu_tab.t <- t(otu_tab)
    seqtab <- matrix(
      as.integer(otu_tab.t),
      nrow=nrow(otu_tab.t),
      ncol=ncol(otu_tab.t),
      dimnames = list(names(derepFs), features)
      )

    # stats of reads
    track <- out
    colnames(track) <- c("input", "filtered")
    rownames(track) <- sample_names
  }

  # Save outputs
  saveRDS(seqtab, file = paste0("rep-seqs_batch_", i, ".rds"))
  saveRDS(track, file = paste0("denoising-stats_batch_", i, ".rds"))
}
if (opt$parallel) {
  # Stops cluster
  parallel::stopCluster(cl)
}

#-------------------------------------------------#
# Collects batches and performs chimera removal   #
#-------------------------------------------------#
denoising_stats <- list.files(path = getwd(),
                              pattern = "denoising-stats_",
                              full.names = TRUE) %>%
  purrr::map(~ data.table::data.table(readRDS(.x))) %>%
  bind_rows()

# Read in seq batches
seqtabs.filenames <- list.files(path = getwd(),
                                pattern = "rep-seqs_",
                                full.names = TRUE)

# Applies mergeSequenceTables on multiple files
if (length(seqtabs.filenames) > 1) {
  seqtabs.list <- lapply(seqtabs.filenames, readRDS)
  seqtabs.merged <- do.call(dada2::mergeSequenceTables, seqtabs.list)
} else {
  seqtabs.merged <- readRDS(seqtabs.filenames)
}

if (!opt$`skip-denoise`) {
  # Remove chimeras
  seqtab.nochim <- dada2::removeBimeraDenovo(
    unqs = seqtabs.merged,
    method = opt$`p-chimera-method`,
    minFoldParentOverAbundance = opt$`p-min-fold-parent-over-abundance`,
    multithread = cpus_n,
    verbose = TRUE
  )
} else {
  seqtab.nochim <- seqtabs.merged
}

#-----------------------------------------#
# Creates dada_report                     #
#-----------------------------------------#

## Adding percentage difference
if (!opt$`skip-denoise`) {
  track <- cbind(denoising_stats, rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoised", "nonchim")
} else {
  track <- denoising_stats
}

track <- add_percentage_columns(track)
track <- cbind(sample = rownames(seqtab.nochim), track)


#-----------------------------------------#
# OUTPUTS FILES                           #
#-----------------------------------------#
# Creating unique fasta file
asv_hash <- unlist(lapply(colnames(seqtab.nochim), function(x) rlang::hash(x)))
dada2::uniquesToFasta(seqtab.nochim, fout="rep-seqs.fna", ids = asv_hash)

# Outputs track
data.table::fwrite(
  x = track,
  file = "denoising-stats.tsv",
  sep = "\t",
  row.names = FALSE
  )

# Creating OTU table
seqtab.nochim <- t(seqtab.nochim) # QIIME has OTUs as rows
col.names <- colnames(seqtab.nochim)
rownames(seqtab.nochim) <- asv_hash
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
utils::write.table(seqtab.nochim, "seq-tab.tsv", sep="\t",
                   row.names=TRUE, col.names=col.names, quote=FALSE)
