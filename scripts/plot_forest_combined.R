#!/usr/bin/env Rscript
# plot_forest_combined.R
# Standalone script to generate a combined forest plot highlighting
# Bentham 2015, Julià 2018, and Meta-Analysis for 13 lead SNPs.

library(data.table)
library(ggplot2)
library(dplyr)

# --- Configuration ---
DIR_RES <- "results"
FILE_META <- file.path(DIR_RES, "sle_meta_lead_snps.tsv")
FILE_LOO  <- file.path(DIR_RES, "leave_one_out_lead_snps.tsv")
FILE_GENE <- file.path(DIR_RES, "gene_annotation.tsv")
OUTPUT_FIG <- file.path(DIR_RES, "forest_plot_combined_discovery.png")

# --- Load Data ---
cat("Loading data...\n")
meta_dt <- fread(FILE_META)
loo_dt  <- fread(FILE_LOO)
gene_dt <- fread(FILE_GENE, header = TRUE, sep = "\t")

# Check column names
cat("Columns in gene_dt:", paste(names(gene_dt), collapse=", "), "\n")

# --- Extract Cohort Estimates from LOO ---
# excl_Julia  -> contains Bentham results
# excl_Bentham -> contains Julià results
cat("Extracting cohort estimates...\n")
bentham <- loo_dt[scenario == "excl_Julia",   .(SNP, BETA = beta_loo, SE = se_loo, Cohort = "Bentham 2015 (UK)")]
julia   <- loo_dt[scenario == "excl_Bentham", .(SNP, BETA = beta_loo, SE = se_loo, Cohort = "Julià 2018 (Spanish)")]
meta    <- meta_dt[, .(SNP, BETA, SE, Cohort = "Meta-Analysis (Combined)")]

# Combine into long format
plot_df <- rbind(bentham, julia, meta)

# Merge gene labels manually to ensure control
cat("Merging gene labels...\n")
if (!"Gene" %in% names(gene_dt)) {
  cat("Warning: 'Gene' column not found in gene_dt. Using SNP as fallback.\n")
  plot_df[, Gene := SNP]
} else {
  plot_df <- merge(plot_df, gene_dt[, .(SNP, Gene)], by="SNP", all.x=TRUE)
}

# Fix empty/NA gene names
plot_df[is.na(Gene) | Gene == "", Gene := SNP]

# Calculate OR and CI
plot_df[, OR := exp(BETA)]
plot_df[, L95 := exp(BETA - 1.96 * SE)]
plot_df[, U95 := exp(BETA + 1.96 * SE)]

# Create Label: Gene (SNP)
plot_df[, Label := paste0(Gene, "\n(", SNP, ")")]

# Order SNPs by Genome Position (CHR and BP) from meta_dt
snp_order <- meta_dt[order(CHR, BP), SNP]
plot_df[, SNP := factor(SNP, levels = rev(snp_order))]

# Get label order
label_order <- unique(plot_df[match(snp_order, SNP), Label])
plot_df[, Label := factor(Label, levels = rev(label_order))]

# Factor for Cohort coloring and layering
plot_df[, Cohort := factor(Cohort, levels = c("Bentham 2015 (UK)", "Julià 2018 (Spanish)", "Meta-Analysis (Combined)"))]

# --- Generate Plot ---
cat("Generating plot...\n")
p <- ggplot(plot_df, aes(x = OR, y = Label, color = Cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.3, alpha = 0.8, 
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(trans = "log", breaks = c(0.5, 0.7, 1, 1.3, 1.5, 2),
                     labels = c("0.5", "0.7", "1.0", "1.3", "1.5", "2.0")) +
  scale_color_manual(values = c("Bentham 2015 (UK)" = "#4393c3", 
                                "Julià 2018 (Spanish)" = "#92c5de", 
                                "Meta-Analysis (Combined)" = "#d6604d")) +
  labs(title = "Combined Forest Plot: Discovery Cohorts vs. Meta-Analysis",
       subtitle = "13 Lead Variants identified in SLE European Meta-Analysis",
       x = "Odds Ratio (OR) and 95% Confidence Interval",
       y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 10)
  )

# Save
cat("Saving to", OUTPUT_FIG, "...\n")
ggsave(OUTPUT_FIG, p, width = 9, height = 11, dpi = 300)

cat("Done.\n")
