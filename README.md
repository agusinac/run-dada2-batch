# Inspiration
Dada2 gets very slow when dealing with large sample size and read depth.
The parallel_dada2.R tries to improve the denoise step when processing large sample size with large read depth on an HPC environment.

This is done by splitting the data into batches and using foreach R module to process samples in parallel with the same learning rates per batch.

The parallel_dada2 works only on single-end reads, since in our organization we merge the paired-end reads with PEAR. We find an enrichment of reads with PEAR versus dada2 denoise paired-end in Qiime2. 

# Mapping file format example (Required):
```
sample-id       absolute-filepath
S103            path/to/S103.fastq.gz
S104            path/to/S104.fastq.gz
```

# Example: Running parallel_dada2.R from command line
```
Rscript --metadata mapping.tsv --batch_n 500 --cpus 8 \
        --p-trunc-q 2 \
        --p-max-ee 6 \
        --p-min-fold-parent-over-abundance 2 \
        --p-chimera-method consensus \
        > dada_report.txt

```
