#!/usr/bin/env Rscript
# scripts/step13_eqtl_plots.R
# Generate multi-tissue dot plot heatmap for SLE eQTL associations.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

DIR_RESULTS <- "results"
FILE_EQTL    <- file.path(DIR_RESULTS, "eqtl_summary.tsv")
OUTPUT_FIG   <- "figures/Fig_4_eQTL_heatmap.png"

cat("=== Generating eQTL Multi-Tissue Visualizations ===\n")

if (!file.exists(FILE_EQTL)) stop("eQTL summary file not found.")

# 1. Load and Clean Data
eqtl <- fread(FILE_EQTL)

# Handle cases where geneSymbol might be missing or ensembl ID
eqtl[, gene_label := ifelse(geneSymbol == "" | is.na(geneSymbol), gencodeId, geneSymbol)]
# Truncate gene names if too long
eqtl[, gene_label := str_trunc(gene_label, 15)]

# Clean Tissue Names
eqtl[, tissue_clean := gsub("_", " ", tissue)]
eqtl[, tissue_clean := gsub("-", " ", tissue_clean)]
eqtl[, tissue_clean := str_to_title(tissue_clean)]

# Filter for readability (Top 30 gene-tissue associations by p-value)
eqtl_top <- eqtl %>%
  arrange(pval) %>%
  head(40)

# 2. Generate Dot Plot Heatmap
cat("Creating dot plot heatmap...\n")

p_eqtl <- ggplot(eqtl_top, aes(x = tissue_clean, y = reorder(gene_label, -pval))) +
  geom_point(aes(size = -log10(pval), fill = nes), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", 
                       midpoint = 0, name = "Normalized Effect Size (NES)\n(Positive = Increased Expression)") +
  scale_size_continuous(range = c(4, 12), name = expression(-log[10](italic(P)[eQTL]))) +
  labs(title = "Multi-Tissue cis-eQTL Landscape of Replicated SLE Loci",
       subtitle = "GTEx v8: Whole Blood, Spleen, and EBV-Lymphocytes",
       x = "Immune-Relevant Tissue", y = "Target Effector Gene") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey30")
  )

# 3. Save Output
ggsave(OUTPUT_FIG, p_eqtl, width = 10, height = 8, dpi = 300, bg = "white")
cat("eQTL visualization complete. Saved to figures/Fig_4_eQTL_heatmap.png\n")
