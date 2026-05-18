#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis"
out_dir <- file.path(base_dir, "2026-05-18_figure_edits")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rep_df <- fread(file.path(base_dir, "results", "spanish_replication_results.tsv"))
master <- fread(file.path(base_dir, "results", "master_results_table.tsv"))

gene_map <- master %>%
  select(RSID, Gene, Effector_Gene) %>%
  mutate(
    display_gene = case_when(
      !is.na(Gene) & Gene != "" ~ Gene,
      !is.na(Effector_Gene) & Effector_Gene != "" ~ Effector_Gene,
      TRUE ~ RSID
    )
  ) %>%
  distinct(RSID, .keep_all = TRUE) %>%
  select(RSID, display_gene)

plot_df <- rep_df %>%
  filter(replicated == TRUE) %>%
  left_join(gene_map, by = "RSID") %>%
  mutate(
    display_gene = ifelse(is.na(display_gene) | display_gene == "", RSID, display_gene),
    order_metric = -log10(P_meta)
  ) %>%
  arrange(desc(order_metric)) %>%
  mutate(
    label = paste0(display_gene, " (", RSID, ")"),
    label = factor(label, levels = rev(unique(label))),
    SE_disc = SE_meta,
    CI_disc_low = BETA_meta - 1.96 * SE_disc,
    CI_disc_high = BETA_meta + 1.96 * SE_disc,
    SE_rep_calc = abs(BETA_rep / qnorm(P_rep / 2, lower.tail = FALSE)),
    CI_rep_low = BETA_rep - 1.96 * SE_rep_calc,
    CI_rep_high = BETA_rep + 1.96 * SE_rep_calc
  )

forest_df <- bind_rows(
  plot_df %>%
    transmute(
      label,
      cohort = "discovery meta-analysis",
      beta = BETA_meta,
      ci_low = CI_disc_low,
      ci_high = CI_disc_high,
      p_value = P_meta
    ),
  plot_df %>%
    transmute(
      label,
      cohort = "Spanish replication",
      beta = BETA_rep,
      ci_low = CI_rep_low,
      ci_high = CI_rep_high,
      p_value = P_rep
    )
) %>%
  mutate(
    cohort = factor(cohort, levels = c("discovery meta-analysis", "Spanish replication")),
    replication_strength = ifelse(p_value < 0.005, "Bonferroni-level", "nominal")
  )

plot_replication_forest <- ggplot(forest_df, aes(x = beta, y = label, color = cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, color = "#7c8798") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    position = position_dodge(width = 0.62),
    width = 0.15,
    linewidth = 0.85,
    orientation = "y"
  ) +
  geom_point(
    aes(shape = replication_strength),
    position = position_dodge(width = 0.62),
    size = 3.05,
    stroke = 0.75
  ) +
  scale_color_manual(
    name = "Dataset",
    values = c(
      "discovery meta-analysis" = "#2c6aa6",
      "Spanish replication" = "#c86418"
    )
  ) +
  scale_shape_manual(
    name = "Replication strength",
    values = c(
      "Bonferroni-level" = 18,
      "nominal" = 16
    )
  ) +
  scale_x_continuous(
    breaks = seq(-0.50, 0.75, by = 0.25),
    limits = c(-0.55, 0.78)
  ) +
  labs(
    title = "Replication forest plot for top reproduced SLE loci",
    x = "Effect size (beta)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0, color = "#1f2937"),
    axis.text.y = element_text(face = "bold", size = 7.6, color = "#1f2937"),
    axis.text.x = element_text(size = 8, color = "#334155"),
    axis.title.x = element_text(size = 8.6, color = "#334155", margin = margin(t = 8)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e5e7eb", linewidth = 0.35),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 0.8, size = 2.8)),
    shape = guide_legend(order = 2)
  )

png_path <- file.path(out_dir, "Replication_forest_top_reproduced_SLE_loci_sentence_case_no_subheading.png")
pdf_path <- file.path(out_dir, "Replication_forest_top_reproduced_SLE_loci_sentence_case_no_subheading.pdf")

ggsave(png_path, plot_replication_forest, width = 7.5, height = 4.8, dpi = 450, bg = "white")
ggsave(pdf_path, plot_replication_forest, width = 7.5, height = 4.8, bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
