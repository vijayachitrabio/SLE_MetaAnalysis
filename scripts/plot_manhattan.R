#!/usr/bin/env Rscript
# plot_manhattan.R
# Standalone script to generate Manhattan plot from meta-analysis results.

library(data.table)
library(qqman)

# --- Configuration ---
INPUT_FILE <- "results/sle_meta_analysis_results.tsv"
OUTPUT_PNG <- "results/sle_meta_manhattan.png"
OUTPUT_PDF <- "results/sle_meta_manhattan.pdf"
TITLE <- "SLE Meta-Analysis (European-Only)"

# Ensure results directory exists
if (!dir.exists("results")) {
  dir.create("results")
}

cat("Loading results from", INPUT_FILE, "...\n")
if (!file.exists(INPUT_FILE)) {
  stop("Input file not found. Please run the meta-analysis script first.")
}

results <- fread(INPUT_FILE)

# Filter for plotting (remove NAs)
plot_dt <- results[!is.na(P) & !is.na(CHR) & !is.na(BP)]
plot_dt[, CHR := as.numeric(CHR)]
plot_dt[, BP := as.numeric(BP)]

cat("Generating Manhattan plot to PNG format...\n")
png(OUTPUT_PNG, width=1200, height=800, res=150)
manhattan(plot_dt, chr="CHR", bp="BP", snp="SNP", p="P", 
          col=c("royalblue", "darkorange"),
          main=TITLE)
dev.off()

cat("Generating Manhattan plot to PDF format...\n")
pdf(OUTPUT_PDF, width=10, height=7)
manhattan(plot_dt, chr="CHR", bp="BP", snp="SNP", p="P", 
          col=c("royalblue", "darkorange"),
          main=TITLE)
dev.off()

cat("Done.\n")
