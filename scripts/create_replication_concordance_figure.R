#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})



rep_df <- fread("results/spanish_replication_results.tsv")

rep_df <- rep_df %>%
  mutate(
    label = case_when(
      !is.na(RSID) & RSID == "rs4853458" ~ "STAT4",
      !is.na(RSID) & RSID == "rs34572943" ~ "ITGAM",
      !is.na(RSID) & RSID == "rs13332649" ~ "IRF8",
      !is.na(RSID) & RSID == "rs10912578" ~ "TNFSF4",
      !is.na(RSID) & RSID == "rs6671847" ~ "FCGR2A",
      !is.na(RSID) & RSID == "rs2647928" ~ "IL12A",
      !is.na(RSID) & RSID == "rs12928726" ~ "CLEC16A",
      TRUE ~ RSID
    ),
    status = case_when(
      replicated ~ "Replicated",
      concordant ~ "Concordant only",
      TRUE ~ "Discordant"
    ),
    discovery_sig = -log10(P_meta),
    replication_sig = -log10(P_rep)
  )

lims <- max(abs(c(rep_df$BETA_meta, rep_df$BETA_rep)), na.rm = TRUE) * 1.15

summary_text <- sprintf(
  "Independent loci: %d\nConcordant effects: %d/%d\nNominal replication: %d/%d\nDiscordant loci: %d",
  nrow(rep_df),
  sum(rep_df$concordant, na.rm = TRUE),
  nrow(rep_df),
  sum(rep_df$replicated, na.rm = TRUE),
  nrow(rep_df),
  sum(!rep_df$concordant, na.rm = TRUE)
)

p <- ggplot(rep_df, aes(x = BETA_meta, y = BETA_rep)) +
  annotate("rect", xmin = 0, xmax = lims, ymin = 0, ymax = lims, fill = "#edf7f0", alpha = 0.7) +
  annotate("rect", xmin = -lims, xmax = 0, ymin = -lims, ymax = 0, fill = "#edf7f0", alpha = 0.7) +
  annotate("rect", xmin = -lims, xmax = 0, ymin = 0, ymax = lims, fill = "#fbefef", alpha = 0.8) +
  annotate("rect", xmin = 0, xmax = lims, ymin = -lims, ymax = 0, fill = "#fbefef", alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "#8b98ab") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "#8b98ab") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7, color = "#6b7280") +
  geom_point(aes(size = discovery_sig, fill = status), shape = 21, color = "#1f2937", stroke = 0.8, alpha = 0.95) +
  geom_text_repel(
    aes(label = label),
    size = 3.7,
    family = "sans",
    color = "#111827",
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    seed = 42,
    max.overlaps = Inf
  ) +
  annotate("label", x = -lims * 0.98, y = lims * 0.98, hjust = 0, vjust = 1,
           label = summary_text, size = 3.7,
           fill = "white", color = "#334155") +
  annotate("text", x = lims * 0.74, y = lims * 0.95, label = "Concordant", color = "#2f6b43", size = 4.0, fontface = "bold") +
  annotate("text", x = lims * 0.72, y = -lims * 0.92, label = "Discordant", color = "#9a3c3c", size = 4.0, fontface = "bold") +
  annotate("text", x = -lims * 0.86, y = -lims * 0.92, label = "Concordant", color = "#2f6b43", size = 4.0, fontface = "bold") +
  scale_fill_manual(values = c(
    "Replicated" = "#2f6db3",
    "Concordant only" = "#7aa6d8",
    "Discordant" = "#d76b6b"
  )) +
  scale_size_continuous(range = c(4, 10), name = expression(-log[10](P[discovery]))) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims), expand = FALSE) +
  labs(
    title = "Discovery-to-Replication Effect Direction Concordance",
    subtitle = "Comparison of discovery meta-analysis and Spanish replication effect sizes across independent SLE loci",
    x = "Discovery effect size (beta)",
    y = "Spanish replication effect size (beta)",
    fill = "Replication status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, color = "#18263d"),
    plot.subtitle = element_text(size = 11, color = "#5b6678"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e5eaf1", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold", color = "#334155"),
    axis.text = element_text(color = "#334155")
  )

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("SLE_Publication_Package/Supplementary/Figures", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/replication_concordance_plot.png", p, width = 9.5, height = 7.5, dpi = 300)
ggsave("figures/replication_concordance_plot.pdf", p, width = 9.5, height = 7.5)
ggsave("SLE_Publication_Package/Supplementary/Figures/replication_concordance_plot.png", p, width = 9.5, height = 7.5, dpi = 300)
ggsave("SLE_Publication_Package/Supplementary/Figures/replication_concordance_plot.pdf", p, width = 9.5, height = 7.5)

message("Saved replication concordance figure to figures/ and supplementary package.")
