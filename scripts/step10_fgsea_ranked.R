#!/usr/bin/env Rscript
# scripts/step10_fgsea_ranked.R
# Advanced Pathway Enrichment using fgsea on Simes-aggregated gene p-values.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(biomaRt)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

# --- Configuration ---
DIR_RESULTS <- "results"
FILE_META    <- file.path(DIR_RESULTS, "discovery_meta_results.tsv")
OUTPUT_NES   <- "figures/pathway_nes_bubbles.png"
OUTPUT_GENE_P <- file.path(DIR_RESULTS, "gene_level_pvalues.tsv")

cat("=== Advanced Pathway Analysis (fgsea) ===\n")

# 1. Fetch Gene Coordinates (GRCh38)
cat("Fetching gene coordinates via biomaRt...\n")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
genes <- getBM(
  attributes = c('external_gene_name', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'chromosome_name',
  values = as.character(1:22),
  mart = ensembl
)
genes <- genes %>% 
  filter(external_gene_name != "") %>%
  mutate(chromosome_name = as.numeric(chromosome_name)) %>%
  rename(CHR = chromosome_name, Start = start_position, End = end_position, Gene = external_gene_name) %>%
  mutate(Start_50kb = Start - 50000, End_50kb = End + 50000)

cat(sprintf("Retrieved %d protein-coding genes.\n", nrow(genes)))

# 2. Load Meta-Analysis Results (Filtered for memory)
cat("Loading meta-analysis results (filtering P < 0.05)...\n")
# We only need SNPs that might contribute to significant gene signals
meta <- fread(FILE_META, select = c("CHR", "BP", "P_meta"))
meta <- meta[P_meta < 0.05]
cat(sprintf("Processing %d SNPs with P < 0.05.\n", nrow(meta)))

# 3. Aggregate SNPs to Genes (Simes method)
cat("Aggregating SNPs to genes (Simes method)...\n")
setDT(genes)
setkey(meta, CHR, BP)

# Rolling join/Overlaps would be better
genes_dt <- genes[, .(Gene, CHR, Start_50kb, End_50kb)]
setnames(genes_dt, c("Start_50kb", "End_50kb"), c("BP_start", "BP_end"))

# Use foverlaps for efficient mapping
meta[, BP_end := BP]
setnames(meta, "BP", "BP_start")
setkey(meta, CHR, BP_start, BP_end)
setkey(genes_dt, CHR, BP_start, BP_end)

overlaps <- foverlaps(genes_dt, meta, nomatch = 0)

# Simes calculation per gene
# P_simes = min( (p_i * n) / i )
gene_pvals <- overlaps[, {
  sorted_p <- sort(P_meta)
  n <- length(sorted_p)
  simes_p <- min( (sorted_p * n) / (1:n) )
  .(SimesP = simes_p)
}, by = Gene]

fwrite(gene_pvals, OUTPUT_GENE_P, sep="\t")
cat("Gene-level p-values saved.\n")

# 4. Run fgsea
cat("Running fgsea with MSigDB Hallmark and Reactome...\n")
# Create ranked list (signed logP if we had direction, but here we use -logP for enrichment)
# Higher rank = more significant SNP evidence for that gene
ranks <- -log10(gene_pvals$SimesP)
names(ranks) <- gene_pvals$Gene
ranks <- sort(ranks, decreasing = TRUE)

# Prepare pathways
msig_h <- msigdbr(species = "human", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
msig_r <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME") %>% split(x = .$gene_symbol, f = .$gs_name)
all_pathways <- c(msig_h, msig_r)

fgsea_res <- fgsea(pathways = all_pathways, stats = ranks, minSize=15, maxSize=500)

# 5. Visualize Top NES
cat("Generating NES bubble plot...\n")
top_pathways <- fgsea_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(20)

# Clean names
top_pathways$pathway <- gsub("HALLMARK_|REACTOME_", "", top_pathways$pathway)
top_pathways$pathway <- gsub("_", " ", top_pathways$pathway)

p_nes <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), color = padj)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Advanced Pathway Enrichment (fgsea)", 
       subtitle = "Simes-aggregated Gene rankings from SLE Discovery Meta-Analysis",
       x = "Normalized Enrichment Score (NES)", y = "")

ggsave(OUTPUT_NES, p_nes, width = 10, height = 7, dpi = 300)
cat("fgsea complete. Plot saved to figures/.\n")
