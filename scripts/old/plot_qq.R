#!/usr/bin/env Rscript
# plot_qq.R
# Standalone script to generate QQ plot and calculate Lambda from meta-analysis results.

library(data.table)
library(qqman)

# --- Configuration ---
INPUT_FILE <- "results/sle_meta_analysis_results.tsv"
OUTPUT_PNG <- "results/sle_meta_qq.png"
OUTPUT_PDF <- "results/sle_meta_qq.pdf"
TITLE_BASE <- "SLE Meta-Analysis: European-only"
NOTE <- "(Bentham GC-corrected pre-meta; Julia uncorrected)"

# Ensure results directory exists
if (!dir.exists("results")) {
  dir.create("results")
}

cat("Loading results from", INPUT_FILE, "...\n")
if (!file.exists(INPUT_FILE)) {
  stop("Input file not found. Please run the meta-analysis script first.")
}

results <- fread(INPUT_FILE)
results <- results[!is.na(P)]

# Calculate Lambda
cat("Calculating Lambda...\n")
chisq <- qchisq(1 - results$P, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
cat("Lambda:", round(lambda, 4), "\n")

# Generate Title
full_title <- paste0(TITLE_BASE, "\nlambda = ", round(lambda, 3), "  ", NOTE)

cat("Generating QQ plot to PNG format...\n")
png(OUTPUT_PNG, width=800, height=800, res=150)
qq(results$P, main=full_title)
dev.off()

cat("Generating QQ plot to PDF format...\n")
pdf(OUTPUT_PDF, width=7, height=7)
qq(results$P, main=full_title)
dev.off()

cat("Done.\n")
