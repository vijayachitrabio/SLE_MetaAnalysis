#!/usr/bin/env Rscript
# scripts/step10_fgsea_ranked.R
# Comprehensive Pathway Enrichment using fgsea on Simes-aggregated gene p-values.
# Databases: Hallmark, Reactome, ImmuneSigDB, GO:BP.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(biomaRt)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

# --- Configuration ---
DIR_RESULTS   <- "results"
DIR_PATHWAY   <- file.path(DIR_RESULTS, "pathway")
if (!dir.exists(DIR_PATHWAY)) dir.create(DIR_PATHWAY, recursive = TRUE)

FILE_META     <- file.path(DIR_RESULTS, "discovery_meta_results.tsv")
OUTPUT_NES    <- "figures/pathway_nes_bubbles.png"
OUTPUT_GENE_P <- file.path(DIR_RESULTS, "gene_level_pvalues.tsv")
OUTPUT_FINAL  <- file.path(DIR_PATHWAY, "fgsea_significant_all.tsv")

cat("=== Advanced Comprehensive Pathway Analysis (fgsea) ===\n")

if (!file.exists(FILE_META)) stop("Meta-analysis results not found.")

# 1. Fetch Gene Coordinates (GRCh38)
cat("Fetching gene coordinates via biomaRt... ")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  attributes = c('external_gene_name', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'chromosome_name',
  values = as.character(1:22),
  mart = ensembl
)

genes <- genes %>% 
  filter(external_gene_name != "") %>%
  filter(!grepl("^HSCHR", chromosome_name)) %>%
  mutate(CHR = as.numeric(chromosome_name)) %>%
  rename(Start = start_position, End = end_position, Gene = external_gene_name) %>%
  mutate(Start_50kb = Start - 50000, End_50kb = End + 50000)

cat(sprintf("Retrieved %d protein-coding genes.\n", nrow(genes)))

# 2. Load Meta-Analysis Results
cat("Loading meta-analysis results (filtering P < 0.05)... ")
meta <- fread(FILE_META, select = c("CHR", "BP", "P_meta"))
meta <- meta[P_meta < 0.05]
cat(sprintf("Processing %d SNPs.\n", nrow(meta)))

# 3. Aggregate SNPs to Genes (Simes method)
cat("Aggregating SNPs to genes (Simes method)... ")
setDT(genes)
setkey(meta, CHR, BP)

genes_dt <- genes[, .(Gene, CHR, Start_50kb, End_50kb)]
setnames(genes_dt, c("Start_50kb", "End_50kb"), c("BP_start", "BP_end"))

meta[, BP_end := BP]
setnames(meta, "BP", "BP_start")
setkey(meta, CHR, BP_start, BP_end)
setkey(genes_dt, CHR, BP_start, BP_end)

overlaps <- foverlaps(genes_dt, meta, nomatch = 0)

gene_pvals <- overlaps[, {
  sorted_p <- sort(P_meta)
  n <- length(sorted_p)
  simes_p <- min( (sorted_p * n) / (1:n) )
  .(SimesP = simes_p)
}, by = Gene]

fwrite(gene_pvals, OUTPUT_GENE_P, sep="\t")
cat("Saved to:", OUTPUT_GENE_P, "\n")

# 4. Run fgsea on multiple databases
cat("Running fgsea across multiple databases...\n")
ranks <- -log10(gene_pvals$SimesP)
names(ranks) <- gene_pvals$Gene
ranks[is.infinite(ranks)] <- max(ranks[is.finite(ranks)]) + 1 # Cap INF
ranks <- sort(ranks, decreasing = TRUE)

# Prepare Database Collections
collections <- list(
  Hallmark = msigdbr(species = "human", category = "H"),
  Reactome = msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME"),
  ImmuneSigDB = msigdbr(species = "human", category = "C7"),
  GO_BP = msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
)

all_results <- list()

for (db_name in names(collections)) {
  cat(sprintf("  Analyzing %s... ", db_name))
  pathways <- split(x = collections[[db_name]]$gene_symbol, f = collections[[db_name]]$gs_name)
  
  res <- fgsea(pathways = pathways, stats = ranks, minSize=15, maxSize=500)
  res$database <- db_name
  
  # Save individual table
  fwrite(res, file.path(DIR_PATHWAY, sprintf("fgsea_%s.tsv", db_name)), sep="\t")
  
  all_results[[db_name]] <- res
  cat("Done.\n")
}

fgsea_all <- rbindlist(all_results)
fwrite(fgsea_all[padj < 0.05], OUTPUT_FINAL, sep="\t")
cat("Consolidated significant results saved.\n")

# 5. Visualize Top NES
cat("Generating NES bubble plot (Top 20 significant)... ")
top_pathways <- fgsea_all %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(20)

# Clean names for plotting
top_pathways$pathway_clean <- gsub("HALLMARK_|REACTOME_|GOBP_", "", top_pathways$pathway)
top_pathways$pathway_clean <- gsub("_", " ", top_pathways$pathway_clean)

p_nes <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_clean, NES), size = -log10(padj), color = padj)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Pathway Enrichment (SLE 47-Loci Architecture)", 
       subtitle = "Ranks: -log10(Simes P-value) from Discovery Meta-Analysis",
       x = "Normalized Enrichment Score (NES)", y = "",
       caption = "Database: Hallmark, Reactome, ImmuneSigDB, GO:BP")

ggsave(OUTPUT_NES, p_nes, width = 11, height = 7, dpi = 300)
cat("Done.\n")
