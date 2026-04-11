#!/usr/bin/env Rscript
# step4_pathway_publication.R
# Publication-quality pathway analysis for SLE European meta-analysis
#
# Strategy (3-layer approach):
#   Layer 1 — MAGMA-style gene-level aggregation (SNP → gene p-value)
#             Uses MAGMA's default window (35 kb upstream, 10 kb downstream)
#             Aggregates all SNPs in a gene window via Simes' method (LD-robust)
#             Output: gene_level_pvalues.tsv
#
#   Layer 2 — Pre-ranked GSEA (fgsea) on gene-level p-values
#             Rank stat = -log10(gene_p) x sign(mean_beta)
#             Gene sets: MSigDB Hallmark + C2 Reactome + C5 GO:BP
#             Output: fgsea_results_*.tsv + forest-style enrichment plots
#
#   Layer 3 — g:Profiler ordered query on full ranked gene list
#             GO:BP, KEGG, Reactome, WikiPathways
#             Correction: g:SCS (more stringent than BH)
#             Output: gprofiler_results.tsv + Manhattan-style plot
#
# NOTE: MAGMA binary not available (server unreachable at time of analysis).
#       Simes' method is the recommended LD-robust alternative for gene-level
#       p-value aggregation without a reference panel requirement. Results are
#       comparable to MAGMA gene-analysis step for pathway input purposes.
#       Reference: Simes RJ (1986) Biometrika; Li & Siegmund (2004) Genetics.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(gprofiler2)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(biomaRt)
})

RESULTS_DIR <- "../results"
FIG_DIR     <- "../figures"
dir.create(file.path(RESULTS_DIR, "pathway"), showWarnings = FALSE)
dir.create(FIG_DIR, showWarnings = FALSE)
PATH_DIR <- file.path(RESULTS_DIR, "pathway")

cat("=== Publication-Quality Pathway Analysis (3-layer) ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# LAYER 1: Gene-level p-value aggregation (MAGMA-style, Simes' method)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Layer 1: Gene-level p-value aggregation (Simes method) ---\n")

meta_dt <- fread(file.path(RESULTS_DIR, "sle_meta_analysis_results.tsv"))
meta_dt <- meta_dt[!is.na(P) & P > 0 & !is.na(CHR) & !is.na(BP)]
meta_dt[, CHR := as.integer(CHR)]
cat("Total SNPs loaded:", nrow(meta_dt), "\n")

# Retrieve gene coordinates from Ensembl GRCh37
cat("Connecting to Ensembl GRCh37 for gene coordinates...\n")
ensembl37 <- tryCatch(
  useMart("ensembl", host = "https://grch37.ensembl.org",
          dataset = "hsapiens_gene_ensembl"),
  error = function(e) { cat("  Ensembl unavailable:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(ensembl37)) {
  gene_coords <- getBM(
    attributes = c("hgnc_symbol", "entrezgene_id",
                   "chromosome_name", "start_position", "end_position", "strand"),
    filters    = "chromosome_name",
    values     = as.character(1:22),
    mart       = ensembl37
  )
  gene_coords <- as.data.table(gene_coords)
  gene_coords <- gene_coords[hgnc_symbol != "" & !is.na(entrezgene_id)]
  gene_coords[, CHR := as.integer(chromosome_name)]
  # MAGMA default window: 35 kb upstream, 10 kb downstream (strand-aware)
  gene_coords[, win_start := ifelse(strand == 1,
                                    start_position - 35000,
                                    start_position - 10000)]
  gene_coords[, win_end   := ifelse(strand == 1,
                                    end_position   + 10000,
                                    end_position   + 35000)]
  gene_coords[win_start < 0, win_start := 0]
  fwrite(gene_coords, file.path(PATH_DIR, "gene_coordinates_grch37.tsv"), sep = "\t")
  cat("Gene coordinates retrieved:", nrow(gene_coords), "genes\n")
} else {
  cat("Using cached gene coordinates if available...\n")
  cache_f <- file.path(PATH_DIR, "gene_coordinates_grch37.tsv")
  if (file.exists(cache_f)) {
    gene_coords <- fread(cache_f)
    cat("Loaded from cache:", nrow(gene_coords), "genes\n")
  } else {
    stop("Ensembl unavailable and no cache found. Run with internet access first.")
  }
}

# Simes' method: for each gene window collect SNP p-values,
# compute gene p = min_k [ k * p_(k) / n ] where p_(k) are ordered p-values
simes_p <- function(pvals) {
  pv <- sort(pvals[!is.na(pvals)])
  n  <- length(pv)
  if (n == 0) return(NA_real_)
  if (n == 1) return(pv[1])
  min((n / seq_along(pv)) * pv)
}

cat("Aggregating SNP p-values to gene level (Simes method)...\n")
setkey(meta_dt, CHR, BP)

gene_results <- gene_coords[, {
  snps_in_window <- meta_dt[CHR == .BY$CHR &
                               BP >= win_start &
                               BP <= win_end]
  if (nrow(snps_in_window) == 0) {
    .(gene_p = NA_real_, n_snps = 0L, mean_beta = NA_real_, best_snp = NA_character_, best_snp_p = NA_real_)
  } else {
    best_idx <- which.min(snps_in_window$P)
    .(gene_p    = simes_p(snps_in_window$P),
      n_snps    = nrow(snps_in_window),
      mean_beta = mean(snps_in_window$BETA, na.rm = TRUE),
      best_snp  = snps_in_window$SNP[best_idx],
      best_snp_p = snps_in_window$P[best_idx])
  }
}, by = .(hgnc_symbol, entrezgene_id, CHR, win_start, win_end)]

gene_results <- gene_results[!is.na(gene_p)]
gene_results[, gene_p_adj := p.adjust(gene_p, method = "BH")]
setorder(gene_results, gene_p)
fwrite(gene_results, file.path(PATH_DIR, "gene_level_pvalues.tsv"), sep = "\t")
cat("Gene-level results:", nrow(gene_results), "genes\n")
cat("Genes with gene_p < 0.05 (BH):", sum(gene_results$gene_p_adj < 0.05, na.rm = TRUE), "\n\n")

# Rank statistic for GSEA: signed -log10(p) captures direction
gene_results[, rank_stat := -log10(gene_p) * sign(mean_beta)]
# Deduplicate: keep highest |rank_stat| per gene symbol (handles multi-mapping)
gene_results_dedup <- gene_results[!is.na(rank_stat) & is.finite(rank_stat)][
  order(-abs(rank_stat))][!duplicated(hgnc_symbol)]
rank_vec <- setNames(gene_results_dedup$rank_stat, gene_results_dedup$hgnc_symbol)
rank_vec <- sort(rank_vec, decreasing = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# LAYER 2: Pre-ranked GSEA with fgsea (MSigDB gene sets)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Layer 2: Pre-ranked GSEA (fgsea) ---\n")

# Load gene sets from MSigDB (msigdbr >= 10.0 uses collection/subcollection)
get_gene_sets <- function(collection, subcollection = NULL) {
  df <- if (is.null(subcollection))
    msigdbr(species = "Homo sapiens", collection = collection)
  else
    msigdbr(species = "Homo sapiens", collection = collection,
            subcollection = subcollection)
  split(df$gene_symbol, df$gs_name)
}

hallmark    <- get_gene_sets("H")
reactome    <- get_gene_sets("C2", "CP:REACTOME")
go_bp       <- get_gene_sets("C5", "GO:BP")
immunologic <- get_gene_sets("C7", "IMMUNESIGDB")

run_fgsea_safe <- function(pathways, ranks, label, min_size = 10, max_size = 500, nperm = 10000) {
  cat("  Running fgsea:", label, "...\n")
  set.seed(42)
  res <- fgsea(pathways  = pathways,
               stats     = ranks,
               minSize   = min_size,
               maxSize   = max_size,
               nPermSimple = nperm)
  res <- as.data.table(res)
  res[, database := label]
  res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
  setorder(res, pval)
  out_f <- file.path(PATH_DIR, paste0("fgsea_", gsub("[^a-zA-Z0-9]", "_", label), ".tsv"))
  fwrite(res, out_f, sep = "\t")
  cat("  Saved:", out_f, "\n")
  cat("  Significant (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
  res
}

fgsea_hallmark    <- run_fgsea_safe(hallmark,    rank_vec, "Hallmark")
fgsea_reactome    <- run_fgsea_safe(reactome,    rank_vec, "Reactome")
fgsea_gobp        <- run_fgsea_safe(go_bp,       rank_vec, "GO_BP")
fgsea_immunologic <- run_fgsea_safe(immunologic, rank_vec, "ImmuneSigDB", max_size = 1000)

# Combine all significant results into one table
all_fgsea <- rbindlist(list(fgsea_hallmark, fgsea_reactome, fgsea_gobp, fgsea_immunologic))
sig_fgsea  <- all_fgsea[padj < 0.05]
setorder(sig_fgsea, pval)
fwrite(sig_fgsea, file.path(PATH_DIR, "fgsea_significant_all.tsv"), sep = "\t")
cat("\nTotal significant pathways across all databases (padj < 0.05):", nrow(sig_fgsea), "\n\n")

# ── fgsea enrichment plot (top 20 per database, publication style) ────────────
plot_fgsea_top <- function(res_dt, label, n = 20) {
  top <- res_dt[padj < 0.05][order(pval)][seq_len(min(.N, n))]
  if (nrow(top) == 0) { cat("  No significant pathways for", label, "— skipping plot\n"); return(invisible(NULL)) }
  top[, pathway_short := stringr::str_wrap(
    gsub("HALLMARK_|REACTOME_|GOBP_|IMMUNESIGDB_", "", pathway), 40)]
  top[, direction := ifelse(NES > 0, "Enriched (SLE risk)", "Depleted (SLE risk)")]

  p <- ggplot(top, aes(x = NES, y = reorder(pathway_short, NES),
                       fill = direction, size = -log10(padj))) +
    geom_point(shape = 21, color = "grey30") +
    scale_fill_manual(values = c("Enriched (SLE risk)" = "#d73027",
                                 "Depleted (SLE risk)"  = "#4575b4")) +
    scale_size_continuous(range = c(3, 10), name = "-log10(adj.P)") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title    = paste0("fgsea: ", label),
         subtitle = paste0("Top ", nrow(top), " significant pathways (FDR < 0.05)"),
         x = "Normalised Enrichment Score (NES)", y = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title    = element_text(face = "bold"),
          axis.text.y   = element_text(size = 8))

  out_f <- file.path(FIG_DIR, paste0("fgsea_", gsub("[^a-zA-Z0-9]", "_", label), "_bubble.png"))
  ggsave(out_f, p, width = 10, height = max(5, nrow(top) * 0.45 + 2), dpi = 200, limitsize = FALSE)
  cat("  Saved plot:", out_f, "\n")
}

cat("--- Generating fgsea enrichment plots ---\n")
plot_fgsea_top(fgsea_hallmark,    "Hallmark")
plot_fgsea_top(fgsea_reactome,    "Reactome")
plot_fgsea_top(fgsea_gobp,        "GO_BP")
plot_fgsea_top(fgsea_immunologic, "ImmuneSigDB")

# ─────────────────────────────────────────────────────────────────────────────
# LAYER 3: g:Profiler ordered query (full ranked gene list)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Layer 3: g:Profiler ordered query ---\n")

# Use ALL genes ranked by gene_p (ascending = most significant first)
ordered_genes <- gene_results[order(gene_p), hgnc_symbol]
ordered_genes <- ordered_genes[!is.na(ordered_genes) & ordered_genes != ""]

cat("Genes submitted to g:Profiler (ordered):", length(ordered_genes), "\n")

gp_res <- tryCatch(
  gost(query            = ordered_genes,
       organism         = "hsapiens",
       ordered_query    = TRUE,          # respects gene ranking by p-value
       significant      = TRUE,
       exclude_iea      = TRUE,          # exclude inferred electronic annotations
       measure_underrepresentation = FALSE,
       evcodes          = TRUE,
       sources          = c("GO:BP", "GO:MF", "KEGG", "REAC", "WP"),
       correction_method = "g_SCS"),     # g:SCS is more stringent than BH
  error = function(e) { cat("g:Profiler error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(gp_res) && !is.null(gp_res$result) && nrow(gp_res$result) > 0) {
  gp_dt <- as.data.table(gp_res$result)
  gp_dt[, intersection := as.character(intersection)]
  fwrite(gp_dt, file.path(PATH_DIR, "gprofiler_results.tsv"), sep = "\t")
  cat("g:Profiler significant terms:", nrow(gp_dt), "\n")

  # Manhattan-style plot (g:Profiler native)
  gp_plot <- gostplot(gp_res, capped = TRUE, interactive = FALSE)
  ggplot2::ggsave(file.path(FIG_DIR, "gprofiler_manhattan.png"),
                  plot = gp_plot, width = 14, height = 7, dpi = 200)
  cat("Saved: figures/gprofiler_manhattan.png\n")

  # Top 20 per source — dotplot
  top_gp <- gp_dt[order(p_value)][, .SD[seq_len(min(.N, 5))], by = source]
  top_gp[, term_name_wrap := stringr::str_wrap(term_name, 45)]
  p_gp <- ggplot(top_gp, aes(x = -log10(p_value),
                               y = reorder(term_name_wrap, -log10(p_value)),
                               color = source, size = intersection_size)) +
    geom_point(alpha = 0.8) +
    scale_color_brewer(palette = "Set1", name = "Database") +
    scale_size_continuous(range = c(3, 10), name = "Genes in overlap") +
    labs(title    = "g:Profiler pathway enrichment (ordered query, g:SCS correction)",
         subtitle = "Top 5 terms per database",
         x = "-log10(adjusted P)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          plot.title    = element_text(face = "bold"),
          axis.text.y   = element_text(size = 9))
  ggsave(file.path(FIG_DIR, "gprofiler_top_dotplot.png"),
         p_gp, width = 12, height = max(6, nrow(top_gp) * 0.4 + 2), dpi = 200)
  cat("Saved: figures/gprofiler_top_dotplot.png\n")

} else {
  cat("No significant g:Profiler results.\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE (methods-ready)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Generating pathway analysis summary ---\n")

summary_rows <- list()

# fgsea summary
for (db in c("Hallmark","Reactome","GO_BP","ImmuneSigDB")) {
  dt <- switch(db,
    "Hallmark"    = fgsea_hallmark,
    "Reactome"    = fgsea_reactome,
    "GO_BP"       = fgsea_gobp,
    "ImmuneSigDB" = fgsea_immunologic
  )
  n_sig <- sum(dt$padj < 0.05, na.rm = TRUE)
  top_path <- if (n_sig > 0) dt[padj < 0.05][order(pval)]$pathway[1] else "—"
  top_nes  <- if (n_sig > 0) round(dt[padj < 0.05][order(pval)]$NES[1], 3) else NA
  summary_rows[[length(summary_rows)+1]] <- data.table(
    Method = "fgsea (pre-ranked GSEA)",
    Database = db,
    N_significant = n_sig,
    Top_pathway = top_path,
    Top_NES = top_nes
  )
}

# g:Profiler summary
if (!is.null(gp_res) && !is.null(gp_res$result) && nrow(gp_res$result) > 0) {
  gp_dt2 <- as.data.table(gp_res$result)
  for (src in unique(gp_dt2$source)) {
    sub <- gp_dt2[source == src][order(p_value)]
    summary_rows[[length(summary_rows)+1]] <- data.table(
      Method = "g:Profiler (ordered query, g:SCS)",
      Database = src,
      N_significant = nrow(sub),
      Top_pathway = sub$term_name[1],
      Top_NES = NA_real_
    )
  }
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(PATH_DIR, "pathway_analysis_summary.tsv"), sep = "\t")
cat("Summary saved to results/pathway/pathway_analysis_summary.tsv\n")
print(summary_dt[, .(Method, Database, N_significant, Top_pathway)])

cat("\n=== Pathway analysis complete ===\n")
cat("Outputs:\n")
cat("  results/pathway/gene_level_pvalues.tsv\n")
cat("  results/pathway/fgsea_*.tsv\n")
cat("  results/pathway/fgsea_significant_all.tsv\n")
cat("  results/pathway/gprofiler_results.tsv\n")
cat("  results/pathway/pathway_analysis_summary.tsv\n")
cat("  figures/fgsea_*_bubble.png\n")
cat("  figures/gprofiler_manhattan.png\n")
cat("  figures/gprofiler_top_dotplot.png\n")
