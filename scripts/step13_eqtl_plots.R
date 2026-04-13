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

# 1. Load and Clean Data
eqtl <- fread(FILE_EQTL)

# Try to map Ensembl IDs that look like Ensembl IDs to gene symbols using biomaRt
if (any(grepl("^ENSG00000", eqtl$geneSymbol))) {
  cat("Attempting to map Ensembl IDs to gene symbols...\n")
  suppressPackageStartupMessages(library(biomaRt))
  tryCatch({
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    ens_ids <- eqtl[grepl("^ENSG00000", geneSymbol), unique(geneSymbol)]
    gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                      filters = "ensembl_gene_id",
                      values = ens_ids,
                      mart = ensembl)
    # Create mapping: use gene name, or readable label for genes without name
    gene_map$geneSymbol_fixed <- ifelse(nchar(gene_map$external_gene_name) > 0,
                                         gene_map$external_gene_name,
                                         ifelse(gene_map$gene_biotype == "TEC", 
                                                "TEC_gene",
                                                paste0("lncRNA_", gene_map$ensembl_gene_id)))
    map_df <- setNames(data.frame(gene_map$ensembl_gene_id, gene_map$geneSymbol_fixed),
                       c("geneSymbol", "geneSymbol_fixed"))
    # Merge and replace
    eqtl <- as.data.table(eqtl)
    eqtl <- merge(eqtl, map_df, by = "geneSymbol", all.x = TRUE)
    eqtl[, geneSymbol := ifelse(is.na(geneSymbol_fixed), geneSymbol, geneSymbol_fixed)]
    eqtl[, geneSymbol_fixed := NULL]
    cat("Mapped", length(ens_ids), "Ensembl IDs\n")
  }, error = function(e) {
    cat("Could not fetch from biomaRt, keeping original IDs\n")
  })
}

# Handle cases where geneSymbol might be missing or ensembl ID
eqtl[, gene_label := ifelse(geneSymbol == "" | is.na(geneSymbol), gencodeId, geneSymbol)]
eqtl[, gene_label := str_trunc(gene_label, 15)]

# Clean Tissue Names
eqtl[, tissue_clean := gsub("_", " ", tissue)]
eqtl[, tissue_clean := gsub("-", " ", tissue_clean)]
eqtl[, tissue_clean := str_to_title(tissue_clean)]

# Add significance indicator
eqtl[, significant := fifelse(pval < 0.05/nrow(eqtl), "FDR < 0.05", ifelse(pval < 0.05, "p < 0.05", "NS"))]

# Filter for readability (Top 40 gene-tissue associations by p-value)
eqtl_top <- eqtl %>%
  arrange(pval) %>%
  head(40)

# 2. Enhanced Dot Plot Heatmap
cat("Creating enhanced dot plot heatmap...\n")

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

ggsave(OUTPUT_DOT, p_eqtl, width = 10, height = 8, dpi = 300, bg = "white")
cat("Dot plot saved to figures/Fig_4_eQTL_heatmap.png\n")

# 3. Clustered Heatmap Matrix (NEW)
cat("Creating clustered heatmap...\n")

# Pivot to matrix format (aggregate duplicates with mean)
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

# Fill missing with 0
nes_matrix[is.na(nes_matrix)] <- 0
pval_matrix[is.na(pval_matrix)] <- 1

# Custom annotation for significance
annotation_col <- data.frame(
  Tissue = colnames(nes_matrix),
  row.names = colnames(nes_matrix)
)
annotation_row <- data.frame(
  Gene = rownames(nes_matrix),
  row.names = rownames(nes_matrix)
)

# Color palette for NES
break_vals <- seq(-1.5, 1.5, length.out = 101)
color_pal <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)

# Create clustered heatmap
pheatmap(nes_matrix,
         color = color_pal,
         breaks = break_vals,
         display_numbers = ifelse(pval_matrix < 0.05, "*", ""),
         number_format = "s",
         fontsize_number = 10,
         cellwidth = 20,
         cellheight = 14,
         treeheight_row = 50,
         treeheight_col = 50,
         main = "Clustered eQTL Heatmap: Tissue-Gene NES",
         legend = TRUE,
         annotation_col = annotation_col,
         filename = OUTPUT_HEAT,
         width = 10,
         height = 8,
         dpi = 300)

cat("Clustered heatmap saved to figures/Fig_4_eQTL_clustered_heatmap.png\n")
cat("eQTL visualization complete.\n")
