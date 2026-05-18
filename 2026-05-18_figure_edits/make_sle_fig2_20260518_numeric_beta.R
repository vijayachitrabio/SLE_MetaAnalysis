#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis"
out_dir <- file.path(base_dir, "2026-05-18_figure_edits")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

top_loci_path <- file.path(base_dir, "results", "top_loci_summary_table.tsv")
master_path <- file.path(base_dir, "results", "master_results_table.tsv")

top_loci_raw <- fread(top_loci_path)
master_results <- fread(master_path)

clean_gene <- function(x) {
  x <- ifelse(is.na(x) | x == "", "Unmapped", x)
  sub("\\s*\\(.*\\)$", "", x)
}

master_gene_labels <- master_results %>%
  mutate(
    LOCUS_ID = paste0(CHR, ":", POS),
    GENE_LABEL = clean_gene(GENE)
  ) %>%
  group_by(LOCUS_ID) %>%
  summarise(GENE_LABEL = first(GENE_LABEL), .groups = "drop")

top_loci <- top_loci_raw %>%
  mutate(
    LOCUS_ID = paste0(CHR, ":", TOP_SNP_POS),
    Locus_class = ifelse(
      NOVEL_LOCUS == "Novel",
      "putative novel",
      "known / nearby known"
    ),
    Support = ifelse(
      LAVA_REPLICATED,
      "LAVA replicated",
      "not replicated"
    ),
    neg_log10_p = -log10(P_MR_RANDOM),
    BETA_ABS = abs(BETA_RANDOM)
  ) %>%
  left_join(master_gene_labels, by = "LOCUS_ID") %>%
  mutate(
    GENE_LABEL = ifelse(is.na(GENE_LABEL), TOP_SNP, GENE_LABEL),
    DISPLAY_LABEL = paste0(GENE_LABEL, "\n", TOP_SNP),
    DISPLAY_LABEL = factor(DISPLAY_LABEL, levels = rev(DISPLAY_LABEL))
  )

plot_fig2 <- ggplot(
  top_loci,
  aes(x = neg_log10_p, y = DISPLAY_LABEL)
) +
  geom_segment(
    aes(x = 0, xend = neg_log10_p, yend = DISPLAY_LABEL, color = Locus_class),
    linewidth = 0.75,
    alpha = 0.75
  ) +
  geom_point(
    aes(size = BETA_ABS, color = Locus_class, shape = Support),
    stroke = 0.55
  ) +
  geom_text(
    aes(label = sprintf("%.2f", BETA_RANDOM)),
    nudge_x = 4.5,
    size = 2.8,
    color = "#3f3f3f",
    hjust = 0
  ) +
  scale_color_manual(
    name = "Locus class",
    values = c(
      "known / nearby known" = "#2f9d78",
      "putative novel" = "#d9601a"
    )
  ) +
  scale_shape_manual(
    name = "Support",
    values = c(
      "LAVA replicated" = 16,
      "not replicated" = 17
    )
  ) +
  scale_size_continuous(name = "beta", range = c(4.2, 8.5)) +
  scale_x_continuous(
    expression(-log[10](italic(P)[meta])),
    limits = c(0, max(top_loci$neg_log10_p) + 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Top 10 SLE loci ranked by statistical strength",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    axis.text.y = element_text(face = "bold", size = 7.7, color = "#2d2d2d"),
    axis.text.x = element_text(size = 8, color = "#4a4a4a"),
    axis.title.x = element_text(size = 8.5, margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e3e3e3", linewidth = 0.35),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3.8, linewidth = 0)),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  )

png_path <- file.path(out_dir, "Fig_2_Top_loci_ranked_numeric_beta_no_subheading.png")
pdf_path <- file.path(out_dir, "Fig_2_Top_loci_ranked_numeric_beta_no_subheading.pdf")

ggsave(png_path, plot_fig2, width = 7.4, height = 4.7, dpi = 450, bg = "white")
ggsave(pdf_path, plot_fig2, width = 7.4, height = 4.7, bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
