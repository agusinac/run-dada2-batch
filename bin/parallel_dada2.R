#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 01-10-2024
# 
#   Optimizes DADA2 in two ways:
#       1. Splits mapping file in batches, can be specified via --batch_n
#       2. Disables multithreading at denoise but allocates each sample on a separate CPU
#       
#   Required INPUT: mapping file
#
#   OUTPUT:
#     1. denoising_stats as RDS file
#     2. seq_tab as RDS file
#
#   Requires Qiime2 to import the data in the right format
#
#   Containers:
#     - parallel_dada2.sif as singularity image
#     - parallel_dada2:latest:cd8a1dd5c5e1 as docker image
#
#------------------------------------------------------------------------------------#

# required libraries
library("foreach")
library("dplyr")

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
                                           help = "sets seed number")
)

# Collects arguments
parser <- optparse::OptionParser(option_list = option_list)
arguments <- optparse::parse_args(parser, positional_arguments=TRUE)

opt <- arguments$options

#-----------------------------------------#
# Loads required files                    #
#-----------------------------------------#

if (!is.null(opt$metadata)) {
  mapping <- data.table::fread(opt$metadata)
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

#-----------------------------------------#
# Setting up parallel and seed            #
#-----------------------------------------#

# Setting up seed
set.seed(seed = seed_n)

# Setting up parallel backend
cl <- parallel::makeCluster(cpus_n)
doParallel::registerDoParallel(cl)

#-----------------------------------------#
# Preparation of batches                  #
#-----------------------------------------#

# Shuffling rows
mapping[base::sample(.N)]
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
  sample_names <- batches[[i]][["SAMPLE-ID"]]
  sample_fastq <- batches[[i]][["absolute-filepath"]]
  
  # Setting filtered paths
  filtFs <- file.path(filtpath, basename(sample_fastq))
  
  # Filtering script
  out <- dada2::filterAndTrim(sample_fastq, filtFs,
                              truncLen = 0,
                              maxEE = 6,
                              truncQ = 2,
                              rm.phix = TRUE,
                              compress = TRUE,
                              verbose = TRUE,
                              multithread = cpus_n,
                              trimLeft = 0)
  
  # Dereplication
  derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
  names(derepFs) <- sample_names
  
  # Learn error rates
  err <- dada2::learnErrors(derepFs, multithread = cpus_n)
  
  # Parallel Denoising
  dds <- foreach::foreach(sam = sample_names, .combine = "c", .packages = "dada2") %dopar% {
    cat("Processing:", sam, "\n")
    list(sam = dada2::dada(derepFs[[ sam ]],
                           err = err,
                           multithread = FALSE))
  }
  # Create sequence table
  dds <- dds[!sapply(dds, is.null)]
  seqtab <- dada2::makeSequenceTable(dds)
  rownames(seqtab) <- sample_names
  
  # stats of reads
  track <- cbind(out, sapply(dds, getN))
  colnames(track) <- c("input", "filtered", "denoised")
  rownames(track) <- sample_names
  
  # Save outputs
  saveRDS(seqtab, file = paste0("rep-seqs_batch_", i, ".rds"))
  saveRDS(track, file = paste0("denoising-stats_batch_", i, ".rds"))
}
# Stops cluster
parallel::stopCluster(cl)

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
 

# Remove chimeras
seqtab.nochim <- dada2::removeBimeraDenovo(seqtabs.merged,
                                           method = "consensus", 
                                           minFoldParentOverAbundance = 2, 
                                           multithread = cpus_n, 
                                           verbose = TRUE)

#-----------------------------------------#
# Creates dada_report                     #
#-----------------------------------------#

track <- cbind(denoising_stats, rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- rownames(seqtab.nochim)

#-----------------------------------------#
# OUTPUTS FILES                           #
#-----------------------------------------#
# Creating unique fasta file
dada2::uniquesToFasta(seqtab.nochim, fout="rep-seqs.fna", ids=colnames(seqtab.nochim))

# Outputs track
utils::write.table(track, file = paste0("denoising-stats.tsv"), sep = "\t")

# Creating OTU table
seqtab.nochim <- t(seqtab.nochim) # QIIME has OTUs as rows
col.names <- colnames(seqtab.nochim)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
utils::write.table(seqtab.nochim, "seq-tab.tsv", sep="\t",
                   row.names=TRUE, col.names=col.names, quote=FALSE)
