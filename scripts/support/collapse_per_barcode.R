#!/usr/bin/env Rscript

# =========================================================
# collapse_per_barcode.R
# Usage:
# Rscript collapse_per_barcode.R <input_dir> <output_file> [min_total_reads]
# - input_dir: folder with per-chromosome TSVs
# - output_file: collapsed output TSV
# - min_total_reads (optional): minimum reads to keep a barcode
# =========================================================

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript collapse_per_barcode.R <input_dir> <output_file> [min_total_reads]")
}

input_dir <- args[1]
output_file <- args[2]
min_total <- ifelse(length(args) >= 3, as.numeric(args[3]), 0)

library(dplyr)
library(readr)

# List all per-chromosome TSVs in the folder
tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
if(length(tsv_files) == 0){
  stop("No TSV files found in input_dir")
}

# Read and combine all TSVs
message("Reading TSVs...")
all_data <- lapply(tsv_files, function(f){
  read_tsv(f, col_types = cols(
    CB = col_character(),
    total_reads = col_double(),
    exonic_reads = col_double(),
    intronic_reads = col_double(),
    intergenic_reads = col_double(),
    mito_reads = col_double()
  ))
})

combined <- bind_rows(all_data)

# Collapse by barcode
message("Collapsing counts per barcode...")
collapsed <- combined %>%
  group_by(CB) %>%
  summarise(
    total_reads = sum(total_reads),
    exonic_reads = sum(exonic_reads),
    intronic_reads = sum(intronic_reads),
    intergenic_reads = sum(intergenic_reads),
    mito_reads = sum(mito_reads),
    .groups = "drop"
  ) %>%
  # Optional: filter low-read barcodes
  filter(total_reads >= min_total) %>%
  # Add derived metrics
  mutate(
    pct_intergenic = intergenic_reads / total_reads,
    pct_intronic = intronic_reads / total_reads,
    pct_mito = mito_reads / total_reads,
    pct_exonic = exonic_reads / total_reads
  )

# Write collapsed table
message("Writing output to ", output_file)
write_tsv(collapsed, output_file)

message("Done!")
