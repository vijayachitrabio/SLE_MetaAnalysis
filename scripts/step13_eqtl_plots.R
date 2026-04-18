#!/usr/bin/env Rscript
# scripts/step13_eqtl_plots.R
# Generate multi-tissue dot plot heatmap for SLE eQTL associations.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(pheatmap)
  library(tidyr)
  library(tibble)
})

DIR_RESULTS <- "results"
FILE_EQTL    <- file.path(DIR_RESULTS, "eqtl_summary.tsv")
OUTPUT_DOT   <- "figures/Fig_4_eQTL_heatmap.png"
OUTPUT_HEAT  <- "figures/Fig_4_eQTL_clustered_heatmap.png"

cat("=== Generating eQTL Multi-Tissue Visualizations ===\n")

if (!file.exists(FILE_EQTL)) stop("eQTL summary file not found.")

## 1. Load and Clean Data
eqtl <- fread(FILE_EQTL)

# Manual Fallback Mapping for common SLE Ensembl IDs if biomaRt is failing
manual_map <- data.frame(
  gcode = c("ENSG00000255310.2", "ENSG00000269918.1", "ENSG00000284957.1", "ENSG00000255310", "ENSG00000269918", "ENSG00000284957"),
  symbol = c("XKR6-AS1", "XKR6", "AL035661.1", "XKR6-AS1", "XKR6", "AL035661.1")
)

# Apply manual map first
eqtl <- merge(eqtl, manual_map, by.x = "geneSymbol", by.y = "gcode", all.x = TRUE)
eqtl[, geneSymbol := ifelse(!is.na(symbol), symbol, geneSymbol)]
eqtl[, symbol := NULL]

# Try to map remaining Ensembl IDs via biomaRt
if (any(grepl("^ENSG00000", eqtl$geneSymbol))) {
  cat("Attempting to map remaining Ensembl IDs to gene symbols...\n")
  suppressPackageStartupMessages(library(biomaRt))
  tryCatch({
    # Use useEnsembl with a specific mirror for better stability
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
    ens_ids <- eqtl[grepl("^ENSG00000", geneSymbol), unique(geneSymbol)]
    gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ens_ids,
                      mart = ensembl)
    
    map_df <- setNames(data.frame(gene_map$ensembl_gene_id, gene_map$external_gene_name),
                       c("geneSymbol", "geneSymbol_fixed"))
    
    eqtl <- merge(as.data.table(eqtl), map_df, by = "geneSymbol", all.x = TRUE)
    eqtl[, geneSymbol := ifelse(is.na(geneSymbol_fixed) | geneSymbol_fixed == "", geneSymbol, geneSymbol_fixed)]
    eqtl[, geneSymbol_fixed := NULL]
  }, error = function(e) {
    cat("Could not fetch from biomaRt mirror, falling back to original IDs\n")
  })
}

# Handle cases where geneSymbol might be missing or ensembl ID
eqtl[, gene_label := ifelse(geneSymbol == "" | is.na(geneSymbol), gencodeId, geneSymbol)]
eqtl[, gene_label := str_trunc(gene_label, 15)]

# Clean Tissue Names
eqtl[, tissue_clean := gsub("_", " ", tissue)]
eqtl[, tissue_clean := gsub("-", " ", tissue_clean)]
eqtl[, tissue_clean := str_to_title(tissue_clean)]

# 2. Enhanced Dot Plot Heatmap
cat("Creating high-quality dot plot heatmap...\n")
eqtl_top <- eqtl[order(pval)][head(1:.N, 40)]

p_eqtl <- ggplot(eqtl_top, aes(x = tissue_clean, y = reorder(gene_label, -pval))) +
  geom_point(aes(size = -log10(pval), fill = nes), shape = 21, color = "black", stroke = 0.6) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", 
                       midpoint = 0, name = "NES") +
  scale_size_continuous(range = c(5, 15), name = expression(-log[10](P))) +
  labs(title = NULL, x = "Tissue", y = "Effector Gene") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.title = element_text(face = "bold", size = 15, color = "black"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3)
  )

ggsave(OUTPUT_DOT, p_eqtl, width = 11, height = 9, dpi = 600, bg = "white")
ggsave(gsub(".png", ".pdf", OUTPUT_DOT), p_eqtl, width = 11, height = 9)

# 3. Clustered Heatmap Matrix
cat("Creating high-quality clustered heatmap...\n")

nes_matrix <- eqtl_top %>%
  dplyr::select(gene_label, tissue_clean, nes) %>%
  dplyr::group_by(gene_label, tissue_clean) %>%
  dplyr::summarise(nes = mean(nes, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tissue_clean, values_from = nes) %>%
  tibble::column_to_rownames("gene_label") %>%
  as.matrix()

pval_matrix <- eqtl_top %>%
  dplyr::select(gene_label, tissue_clean, pval) %>%
  dplyr::group_by(gene_label, tissue_clean) %>%
  dplyr::summarise(pval = min(pval, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tissue_clean, values_from = pval) %>%
  tibble::column_to_rownames("gene_label") %>%
  as.matrix()

nes_matrix[is.na(nes_matrix)] <- 0
pval_matrix[is.na(pval_matrix)] <- 1

# Pheatmap Generation
pheatmap(nes_matrix,
         color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
         breaks = seq(-1.5, 1.5, length.out = 101),
         display_numbers = ifelse(pval_matrix < 0.05/nrow(eqtl), "**", ifelse(pval_matrix < 0.05, "*", "")),
         fontsize_number = 14,
         number_color = "black",
         cellwidth = 25,
         cellheight = 18,
         fontsize_row = 12,
         fontsize_col = 12,
         angle_col = 45,
         main = "", # Remove heading
         legend = TRUE,
         filename = OUTPUT_HEAT,
         width = 11,
         height = 9)

# Save PDF version
pheatmap(nes_matrix,
         color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
         breaks = seq(-1.5, 1.5, length.out = 101),
         display_numbers = ifelse(pval_matrix < 0.05/nrow(eqtl), "**", ifelse(pval_matrix < 0.05, "*", "")),
         fontsize_number = 14,
         number_color = "black",
         cellwidth = 25,
         cellheight = 18,
         fontsize_row = 12,
         fontsize_col = 12,
         angle_col = 45,
         main = "",
         legend = TRUE,
         filename = gsub(".png", ".pdf", OUTPUT_HEAT),
         width = 11,
         height = 9)

cat("Clustered heatmap saved to figures/Fig_4_eQTL_clustered_heatmap.png\n")
cat("eQTL visualization complete.\n")
