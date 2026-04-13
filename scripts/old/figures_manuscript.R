#!/usr/bin/env Rscript
# figures_manuscript.R
# Generates all publication-ready figures for the SLE European meta-analysis manuscript.
#
# Figures produced:
#   Forest plots
#     Fig_S1_forest_combined.png        — all 13 lead variants, meta effect, replication status
#     Fig_S2_forest_<snp>.png           — per-SNP 4-cohort panels (Bentham / Julià / Meta / FinnGen)
#   Pathway
#     Fig_3_pathway_fgsea_combined.png  — top fgsea pathways across all databases (main figure)
#     Fig_S3_pathway_gprofiler.png      — g:Profiler Manhattan (supplementary)
#     Fig_S4_pathway_gprofiler_dot.png  — g:Profiler top-terms dotplot (supplementary)
#
# All figures saved at 300 DPI, suitable for direct journal submission.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(gprofiler2)
  library(stringr)
})

RES  <- "../results"
FIG  <- "../figures"
dir.create(FIG, showWarnings = FALSE)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "grey40"),
      axis.text     = element_text(color = "black"),
      axis.title    = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title  = element_text(face = "bold", size = base_size - 1),
      legend.text   = element_text(size = base_size - 1),
      strip.text    = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
}

pal_rep <- c("Replicated" = "#2166ac", "Not replicated" = "#d73027",
             "Palindromic\nexcluded" = "grey60", "Not tested" = "grey80")

# ─────────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
cat("Loading data...\n")

lead   <- fread(file.path(RES, "sle_meta_lead_snps.tsv"))
annot  <- fread(file.path(RES, "gene_annotation.tsv"),
                select = c("SNP","Gene","Region","Known_SLE"))
annot[is.na(Gene) | Gene == "", Gene := "Unannotated"]

het    <- fread(file.path(RES, "heterogeneity_lead_snps.tsv"))
fg     <- fread(file.path(RES, "finngen_replication_v2_results_full.tsv"))
loo    <- fread(file.path(RES, "leave_one_out_lead_snps.tsv"))
conc   <- fread(file.path(RES, "european_replication_concordance_table.tsv"))

# Reconstruct per-cohort effects from LOO
# excl_Julia  → remaining = Bentham
# excl_Bentham → remaining = Julià
bentham <- loo[scenario == "excl_Julia",   .(SNP, BETA_b = beta_loo, SE_b = se_loo)]
julia   <- loo[scenario == "excl_Bentham", .(SNP, BETA_j = beta_loo, SE_j = se_loo)]

# FinnGen (exclude palindromic)
fg_sub <- fg[, .(SNP, BETA_fg, SE_fg, P_fg, align_status, replicated)]

# Master table
master <- lead[, .(SNP, CHR, BP, A1, A2, BETA_m = BETA, SE_m = SE, P_m = P)]
master <- merge(master, annot,   by = "SNP", all.x = TRUE)
master <- merge(master, het[, .(SNP, I2, p_Q)], by = "SNP", all.x = TRUE)
master <- merge(master, bentham, by = "SNP", all.x = TRUE)
master <- merge(master, julia,   by = "SNP", all.x = TRUE)
master <- merge(master, fg_sub,  by = "SNP", all.x = TRUE)
master <- merge(master, conc[, .(SNP, replication_status)], by = "SNP", all.x = TRUE)

master[, rep_label := fcase(
  replication_status == "replicated",                       "Replicated",
  replication_status == "palindromic_excluded",             "Palindromic\nexcluded",
  grepl("not_replicated", replication_status),              "Not replicated",
  default = "Not tested"
)]

# Label: Gene (chromosome:position)
master[, label := ifelse(Gene == "Unannotated",
                         paste0("Chr", CHR, ":", format(BP, big.mark=",")),
                         Gene)]
master[, label_full := paste0(label, "\n", A1, ">", A2)]

setorder(master, CHR, BP)
master[, snp_order := .I]   # genomic order index

cat("Lead variants loaded:", nrow(master), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Combined forest plot — all 13 variants, meta effect
# ─────────────────────────────────────────────────────────────────────────────
cat("Generating Fig_1_forest_combined.png ...\n")

# CI
master[, lo_m := BETA_m - 1.96 * SE_m]
master[, hi_m := BETA_m + 1.96 * SE_m]

# Annotation text for right margin
master[, annot_text := paste0(
  "OR=", round(exp(BETA_m), 2),
  " (", round(exp(lo_m), 2), "–", round(exp(hi_m), 2), ")",
  "  P=", formatC(P_m, format = "e", digits = 1),
  "  I²=", I2, "%"
)]

fig1 <- ggplot(master, aes(y = reorder(label_full, -snp_order))) +
  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # CIs
  geom_errorbarh(aes(xmin = exp(lo_m), xmax = exp(hi_m), color = rep_label),
                 height = 0.25, linewidth = 0.7) +
  # Points (size = -log10 P)
  geom_point(aes(x = exp(BETA_m), color = rep_label, size = -log10(P_m)),
             shape = 18) +
  # Replication symbol on right
  geom_text(aes(x = max(exp(hi_m), na.rm = TRUE) * 1.05,
                label = ifelse(rep_label == "Replicated", "\u2713",
                        ifelse(rep_label == "Palindromic\nexcluded", "P",
                        ifelse(rep_label == "Not replicated", "\u2717", ""))),
                color = rep_label),
            size = 4, hjust = 0, fontface = "bold") +
  scale_color_manual(values = pal_rep, name = "FinnGen R12") +
  scale_size_continuous(range = c(3, 8), name = expression(-log[10](P[meta]))) +
  scale_x_continuous(name = "Odds Ratio (95% CI)",
                     trans = scales::log_trans(),
                     breaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 1.8),
                     labels = c("0.6","0.7","0.8","0.9","1.0","1.1","1.3","1.5","1.8"),
                     limits = c(0.55, 2.0),
                     oob    = scales::squish) +
  labs(title    = "SLE European Meta-Analysis: 13 Lead Variants",
       subtitle  = paste0("Discovery: Bentham 2015 + Julià 2018 (N=4,835 cases / 8,517 controls)",
                          "\nReplication: FinnGen R12 L12_LUPUS (~850 cases, Finnish)"),
       y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y  = element_text(size = 9, face = "italic"),
        legend.box   = "vertical")

ggsave(file.path(FIG, "Fig_1_forest_combined.png"),
       fig1, width = 10, height = 7, dpi = 300)
cat("  Saved: figures/Fig_1_forest_combined.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE S1: Per-SNP 4-cohort forest plots (all 13, multi-panel)
# ─────────────────────────────────────────────────────────────────────────────
cat("Generating Fig_S1_forest_per_snp.png (multi-panel) ...\n")

plot_one_snp <- function(row) {
  snp_id <- row$SNP

  plot_df <- data.frame(
    cohort = c("Bentham 2015\n(UK European)",
               "Julià 2018\n(Spanish European)",
               "Meta-analysis\n(IVW fixed-effects)",
               "FinnGen R12\n(Finnish, replication)"),
    BETA   = c(row$BETA_b, row$BETA_j, row$BETA_m, row$BETA_fg),
    SE     = c(row$SE_b,   row$SE_j,   row$SE_m,   row$SE_fg),
    type   = c("Discovery","Discovery","Meta","Replication")
  ) %>%
    filter(!is.na(BETA), !is.na(SE)) %>%
    mutate(
      OR    = exp(BETA),
      lo    = exp(BETA - 1.96 * SE),
      hi    = exp(BETA + 1.96 * SE),
      cohort = factor(cohort, levels = rev(cohort))
    )

  gene_label <- if (row$Gene == "Unannotated")
    paste0("Chr", row$CHR, ":", format(row$BP, big.mark=","))
  else row$Gene

  subtitle_txt <- paste0(
    row$A1, ">", row$A2,
    "  |  Meta OR=", round(exp(row$BETA_m), 2),
    " (P=", formatC(row$P_m, format="e", digits=1), ")",
    "  |  I²=", row$I2, "%"
  )

  ggplot(plot_df, aes(x = OR, y = cohort, color = type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.7) +
    geom_point(aes(size = type), shape = 18) +
    scale_color_manual(values = c("Discovery"   = "#4393c3",
                                  "Meta"        = "#d6604d",
                                  "Replication" = "#74c476"),
                       guide = "none") +
    scale_size_manual(values = c("Discovery" = 3.5, "Meta" = 5, "Replication" = 3.5),
                      guide = "none") +
    scale_x_continuous(trans = "log",
                       name = "Odds Ratio (95% CI)") +
    labs(title    = paste0(snp_id, "  [", gene_label, "]"),
         subtitle = subtitle_txt,
         y = NULL) +
    theme_pub(base_size = 9) +
    theme(plot.title    = element_text(size = 9,  face = "bold.italic"),
          plot.subtitle = element_text(size = 7.5, color = "grey30"),
          axis.text.y   = element_text(size = 8))
}

plots <- lapply(seq_len(nrow(master)), function(i) plot_one_snp(master[i]))

# 3-column grid
combined_panel <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title    = "Supplementary Fig S1 — Per-variant forest plots",
    subtitle = "Cohorts: Bentham 2015 (UK), Julià 2018 (Spanish), Meta-analysis, FinnGen R12 (replication)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(FIG, "Fig_S1_forest_per_snp.png"),
       combined_panel, width = 18, height = 20, dpi = 300, limitsize = FALSE)
cat("  Saved: figures/Fig_S1_forest_per_snp.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2: fgsea — top significant pathways, combined across databases
# ─────────────────────────────────────────────────────────────────────────────
cat("Generating Fig_2_pathway_fgsea.png ...\n")

fgsea_all <- fread(file.path(RES, "pathway/fgsea_significant_all.tsv"))

# Clean pathway names
fgsea_all[, pathway_clean := pathway]
fgsea_all[, pathway_clean := gsub("^HALLMARK_",   "",  pathway_clean)]
fgsea_all[, pathway_clean := gsub("^REACTOME_",   "",  pathway_clean)]
fgsea_all[, pathway_clean := gsub("^GOBP_",       "GO: ", pathway_clean)]
fgsea_all[, pathway_clean := gsub("^GSE[^_]+_",   "",  pathway_clean)]
fgsea_all[, pathway_clean := gsub("_", " ", pathway_clean)]
fgsea_all[, pathway_clean := str_to_sentence(pathway_clean)]
fgsea_all[, pathway_clean := str_wrap(pathway_clean, 42)]

# Top N per database (avoid overcrowding)
top_n_per_db <- fgsea_all[padj < 0.05][order(pval)][
  , .SD[seq_len(min(.N, 6))], by = database]

top_n_per_db[, db_label := factor(database,
  levels = c("Hallmark","Reactome","GO_BP","ImmuneSigDB"),
  labels = c("MSigDB Hallmark","Reactome","GO: Biological Process","ImmuneSigDB"))]

top_n_per_db[, direction := ifelse(NES > 0, "Risk-increasing", "Risk-decreasing")]

fig2 <- ggplot(top_n_per_db,
               aes(x = NES,
                   y = reorder(pathway_clean, NES),
                   fill = direction,
                   size = -log10(padj))) +
  geom_point(shape = 21, color = "grey25", alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  facet_wrap(~ db_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Risk-increasing"  = "#b2182b",
                               "Risk-decreasing"  = "#2166ac"),
                    name = "Direction of effect") +
  scale_size_continuous(range = c(3, 10),
                        name  = expression(-log[10](adj.~italic(P)))) +
  labs(title    = "Gene-set enrichment analysis: SLE-associated pathways",
       subtitle  = paste0("Pre-ranked fgsea on gene-level p-values (Simes method, MAGMA-style window)\n",
                          "Showing top 6 per database (FDR < 0.05)"),
       x = "Normalised Enrichment Score (NES)",
       y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y   = element_text(size = 8.5),
        strip.background = element_rect(fill = "grey92", color = NA),
        legend.box    = "vertical",
        panel.spacing = unit(0.8, "lines"))

ggsave(file.path(FIG, "Fig_2_pathway_fgsea.png"),
       fig2, width = 13, height = 9, dpi = 300)
cat("  Saved: figures/Fig_2_pathway_fgsea.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE S2: g:Profiler — top terms per source, dotplot
# ─────────────────────────────────────────────────────────────────────────────
cat("Generating Fig_S2_pathway_gprofiler_dot.png ...\n")

gp <- fread(file.path(RES, "pathway/gprofiler_results.tsv"))
gp[, term_wrap := str_wrap(term_name, 45)]

# Top 6 per source, order by p-value
top_gp <- gp[order(p_value)][, .SD[seq_len(min(.N, 6))], by = source]
top_gp[, source_label := factor(source,
  levels = c("GO:BP","GO:MF","KEGG","REAC","WP"),
  labels = c("GO: Biological Process","GO: Molecular Function",
             "KEGG","Reactome","WikiPathways"))]

db_colors <- c("GO: Biological Process" = "#e41a1c",
               "GO: Molecular Function" = "#ff7f00",
               "KEGG"                   = "#4daf4a",
               "Reactome"               = "#377eb8",
               "WikiPathways"           = "#984ea3")

fig_s2 <- ggplot(top_gp,
                 aes(x = -log10(p_value),
                     y = reorder(term_wrap, -log10(p_value)),
                     color = source_label,
                     size  = intersection_size)) +
  geom_point(alpha = 0.85) +
  geom_segment(aes(xend = 0, yend = term_wrap),
               linewidth = 0.3, alpha = 0.3) +
  scale_color_manual(values = db_colors, name = "Database") +
  scale_size_continuous(range = c(3, 10), name = "Genes in overlap") +
  labs(title    = "Supplementary Fig S2 — g:Profiler pathway enrichment",
       subtitle  = paste0("Ordered query (all 21,597 genes ranked by gene-level p-value)\n",
                          "g:SCS multiple testing correction | IEA annotations excluded | Top 6 per database"),
       x = expression(-log[10](adjusted~italic(P))),
       y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(size = 8.5),
        legend.box  = "vertical")

ggsave(file.path(FIG, "Fig_S2_pathway_gprofiler_dot.png"),
       fig_s2, width = 13, height = 9, dpi = 300)
cat("  Saved: figures/Fig_S2_pathway_gprofiler_dot.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE S3: g:Profiler Manhattan (native plot, high-res)
# ─────────────────────────────────────────────────────────────────────────────
cat("Generating Fig_S3_pathway_gprofiler_manhattan.png ...\n")

# Re-run gostplot from saved results using token (if online) or reconstruct
# Since we have the result data, rebuild the gostplot object from saved TSV
gp_terms <- gp$term_id
gp_pvals <- gp$p_value

# Use publish_gosttable for a clean table figure of top hits
top_gp_pub <- top_gp[order(p_value)][seq_len(min(nrow(top_gp), 30))]

fig_s3 <- ggplot(top_gp_pub,
                 aes(x = reorder(term_wrap, -log10(p_value)),
                     y = -log10(p_value),
                     fill = source_label)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = db_colors, name = "Database") +
  labs(title    = "Supplementary Fig S3 — g:Profiler: significant pathway terms",
       subtitle  = "Dashed line = FDR 0.05 threshold (g:SCS correction)",
       x = NULL,
       y = expression(-log[10](adjusted~italic(P)))) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "right")

ggsave(file.path(FIG, "Fig_S3_pathway_gprofiler_bar.png"),
       fig_s3, width = 13, height = 10, dpi = 300)
cat("  Saved: figures/Fig_S3_pathway_gprofiler_bar.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== All manuscript figures generated ===\n\n")
cat("Main figures (manuscript body):\n")
cat("  Fig_1_forest_combined.png     — combined forest plot, 13 lead variants\n")
cat("  Fig_2_pathway_fgsea.png       — fgsea NES bubble, 4 databases\n\n")
cat("Supplementary figures:\n")
cat("  Fig_S1_forest_per_snp.png     — individual 4-cohort forest panels\n")
cat("  Fig_S2_pathway_gprofiler_dot.png — g:Profiler dotplot, top 6 per source\n")
cat("  Fig_S3_pathway_gprofiler_bar.png — g:Profiler bar chart, all sig. terms\n")
