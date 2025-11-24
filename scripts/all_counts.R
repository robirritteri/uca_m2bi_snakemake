#!/usr/bin/env Rscript

###############################################
# Merge featureCounts outputs
###############################################

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: all_counts.R <counts_dir> <output_file>")
}

counts_dir <- args[1]
out_file   <- args[2]

# List all featureCounts files produced by Snakemake
files <- list.files(counts_dir, pattern = "_counts\\.txt$", full.names = TRUE)

if (length(files) == 0) {
  stop("No *_counts.txt files found in ", counts_dir)
}

read_fc <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Extract sample name from filename: <sample>_counts.txt -> <sample>
  sample <- basename(path)
  sample <- sub("\\.txt$", "", sample)
  sample <- sub("_counts$", "", sample)
  
  # featureCounts format: first column "Geneid", 7th column = count
  df_small <- df[, c("Geneid", colnames(df)[7])]
  colnames(df_small)[2] <- sample
  return(df_small)
}

# Read all files
list_fc <- lapply(files, read_fc)

# Merge by Geneid
merged <- Reduce(function(x, y) merge(x, y, by = "Geneid"), list_fc)

# Write merged table
write.table(
  merged,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
