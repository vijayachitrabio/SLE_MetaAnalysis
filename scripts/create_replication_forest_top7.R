#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})



rep_df <- fread("results/spanish_replication_results.tsv")
master <- fread("results/master_results_table.tsv")

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
      Cohort = "Discovery meta-analysis",
      Beta = BETA_meta,
      CI_low = CI_disc_low,
      CI_high = CI_disc_high,
      P_value = P_meta
    ),
  plot_df %>%
    transmute(
      label,
      Cohort = "Spanish replication",
      Beta = BETA_rep,
      CI_low = CI_rep_low,
      CI_high = CI_rep_high,
      P_value = P_rep
    )
) %>%
  mutate(
    Cohort = factor(Cohort, levels = c("Discovery meta-analysis", "Spanish replication")),
    signif_label = ifelse(P_value < 0.005, "Bonferroni-level", "Nominal")
  )

p <- ggplot(forest_df, aes(x = Beta, y = label, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "#7c8798") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), position = position_dodge(width = 0.62), height = 0.16, linewidth = 0.9) +
  geom_point(aes(shape = signif_label), position = position_dodge(width = 0.62), size = 3.3, stroke = 0.8) +
  scale_color_manual(values = c(
    "Discovery meta-analysis" = "#2b6cb0",
    "Spanish replication" = "#d97706"
  )) +
  scale_shape_manual(values = c(
    "Bonferroni-level" = 18,
    "Nominal" = 16
  )) +
  labs(
    title = "Replication Forest Plot for Top Reproduced SLE Loci",
    subtitle = "Effect sizes and 95% confidence intervals in discovery meta-analysis and Spanish replication",
    x = "Effect size (beta)",
    y = NULL,
    color = "Dataset",
    shape = "Replication strength"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, color = "#18263d"),
    plot.subtitle = element_text(size = 11, color = "#5b6678"),
    axis.text.y = element_text(face = "bold", color = "#1f2937"),
    axis.text.x = element_text(color = "#334155"),
    axis.title.x = element_text(face = "bold", color = "#334155"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("SLE_Publication_Package/Supplementary/Figures", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/replication_forest_top7.png", p, width = 10.5, height = 6.8, dpi = 300)
ggsave("figures/replication_forest_top7.pdf", p, width = 10.5, height = 6.8)
ggsave("SLE_Publication_Package/Supplementary/Figures/replication_forest_top7.png", p, width = 10.5, height = 6.8, dpi = 300)
ggsave("SLE_Publication_Package/Supplementary/Figures/replication_forest_top7.pdf", p, width = 10.5, height = 6.8)

message("Saved replication forest plot to figures/ and supplementary package.")
