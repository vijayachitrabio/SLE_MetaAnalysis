#!/usr/bin/env Rscript
# scripts/generate_publication_summary_stats.R
# Generates clean, publication-ready full summary statistics files (gzipped)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("=== Generating Publication-Ready Summary Statistics ===\n")

# Set working directory
setwd(".")

# 1. Load full meta-analysis results
cat("Loading full meta-analysis results (1.3 GB)...\n")
meta_path <- "results/discovery_meta_results.tsv"
if (!file.exists(meta_path)) {
  stop("Error: discovery_meta_results.tsv not found.")
}

# Read only necessary columns to save memory
meta <- fread(
  meta_path,
  select = c("RSID", "CHR", "BP", "OA", "EA", "BETA_meta", "SE_meta", "P_meta", "I2", "HetP")
)

raw_n <- nrow(meta)
missing_rsid_n <- meta[is.na(RSID) | RSID == "", .N]
missing_required_n <- meta[
  is.na(RSID) | RSID == "" |
    is.na(CHR) | is.na(BP) |
    is.na(OA) | OA == "" |
    is.na(EA) | EA == "" |
    is.na(BETA_meta) | is.na(SE_meta) | is.na(P_meta),
  .N
]

# 2. Format for Publication / GWAS Catalog
cat("Filtering and formatting columns...\n")
pub_stats <- meta[
  !(is.na(RSID) | RSID == "" |
      is.na(CHR) | is.na(BP) |
      is.na(OA) | OA == "" |
      is.na(EA) | EA == "" |
      is.na(BETA_meta) | is.na(SE_meta) | is.na(P_meta))
] %>%
  rename(
    rsid = RSID,
    chromosome = CHR,
    base_pair_location = BP,
    other_allele = OA,
    effect_allele = EA,
    beta = BETA_meta,
    standard_error = SE_meta,
    p_value = P_meta,
    het_i2 = I2,
    het_p_value = HetP
  )

# 3. Create Output Directories
cat("Creating output directories...\n")
target_dirs <- c(
  "SLE_Publication_Package/Supplementary/Data",
  "2026-04-25_chat_gpt_submission_ready/Supplementary/Data"
)

for (d in target_dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# 4. Save Gzipped Results and QC Report
output_filename <- "SLE_MetaAnalysis_Summary_Statistics_Final_Filtered.tsv.gz"
qc_filename <- "SLE_MetaAnalysis_Summary_Statistics_Final_Filtered_QC.tsv"

qc_report <- data.table(
  Metric = c(
    "Source file",
    "Raw variants in source",
    "Variants excluded for missing required fields",
    "Variants excluded for missing RSID",
    "Final variants retained",
    "Genome-wide significant variants (P < 5e-8)",
    "Suggestive variants (P < 1e-5)",
    "Minimum p-value",
    "Median heterogeneity I2",
    "Mean heterogeneity I2"
  ),
  Value = c(
    meta_path,
    raw_n,
    missing_required_n,
    missing_rsid_n,
    nrow(pub_stats),
    sum(pub_stats$p_value < 5e-8, na.rm = TRUE),
    sum(pub_stats$p_value < 1e-5, na.rm = TRUE),
    format(min(pub_stats$p_value, na.rm = TRUE), scientific = TRUE),
    median(pub_stats$het_i2, na.rm = TRUE),
    mean(pub_stats$het_i2, na.rm = TRUE)
  )
)

cat("Writing and compressing full summary statistics (this may take a minute)...\n")
for (d in target_dirs) {
  out_path <- file.path(d, output_filename)
  qc_path <- file.path(d, qc_filename)
  cat(paste("Saving to:", out_path, "...\n"))
  fwrite(pub_stats, out_path, sep = "\t", compress = "gzip")
  cat(paste("Saving QC report to:", qc_path, "...\n"))
  fwrite(qc_report, qc_path, sep = "\t")
}

cat("\n=== Summary Statistics Generation Complete ===\n")
cat("Generated file: ", output_filename, "\n")
cat("Total variants: ", nrow(pub_stats), "\n")
cat("Genome-wide significant variants: ", sum(pub_stats$p_value < 5e-8, na.rm = TRUE), "\n")
cat("Finished.\n")
