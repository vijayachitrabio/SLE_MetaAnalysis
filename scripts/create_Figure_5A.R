#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dt <- fread("results/Supplementary_Table_6.tsv")
setnames(dt, c("term", "source", "gene_count", "p_value", "fdr"))

wanted <- data.table(
  display_order = 1:12,
  term = c(
    "immune system process",
    "immune response",
    "defense response",
    "interleukin-12-mediated signaling pathway",
    "regulation of immune response",
    "cell activation",
    "Immune System",
    "Interleukin-35 Signalling",
    "Cytokine Signaling in Immune system",
    "Interleukin-12 signaling",
    "Interleukin-12 family signaling",
    "Interferon alpha/beta signaling"
  ),
  term_label = c(
    "Immune system process",
    "Immune response",
    "Defense response",
    "Interleukin-12-mediated\nsignaling",
    "Regulation of immune response",
    "Cell activation",
    "Immune system",
    "Interleukin-35 signaling",
    "Cytokine signaling in immune system",
    "Interleukin-12 signaling",
    "Interleukin-12 family signaling",
    "Interferon alpha/beta signaling"
  ),
  source = c(rep("GO:BP", 6), rep("REAC", 6))
)

plot_dt <- merge(wanted, dt, by = c("term", "source"), sort = FALSE)
setorder(plot_dt, display_order)
plot_dt[, source_label := fifelse(source == "GO:BP", "GO biological process", "Reactome")]
plot_dt[, log_p := -log10(fdr)]
plot_dt[, source_label := factor(source_label, levels = c("GO biological process", "Reactome"))]
plot_dt[, term_plot := term_label]
plot_dt[, term_plot := factor(term_plot, levels = rev(term_plot))]

source_colors <- c(
  "GO biological process" = "#2f837a",
  "Reactome" = "#cf741c"
)

p <- ggplot(plot_dt, aes(x = log_p, y = term_plot, fill = source_label)) +
  geom_col(width = 0.78) +
  geom_text(
    aes(label = gene_count),
    hjust = -0.16,
    size = 3.8,
    color = "#333333"
  ) +
  scale_fill_manual(values = source_colors, name = NULL) +
  scale_x_continuous(
    limits = c(0, 2.62),
    breaks = seq(0, 2.5, by = 0.5),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    caption = "Numbers at bar ends indicate gene count."
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.key.size = unit(0.16, "in"),
    legend.text = element_text(size = 9.4, color = "#111111", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e6eaef", linewidth = 0.35),
    axis.text.y = element_text(size = 8.8, color = "#111111", face = "bold"),
    axis.text.x = element_text(size = 8.8, color = "#4b5563", face = "bold"),
    axis.title.x = element_text(size = 10, color = "#111111", face = "bold", margin = margin(t = 6)),
    plot.caption = element_text(size = 8.2, color = "#4b5563", hjust = 0),
    plot.margin = margin(5, 16, 5, 5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_5A.png", p, width = 7.2, height = 4.05, dpi = 600, bg = "white")
ggsave("figures/Figure_5A.tiff", p, width = 7.2, height = 4.05, dpi = 600, bg = "white", compression = "lzw")
ggsave("figures/5A_new.png", p, width = 7.2, height = 4.05, dpi = 600, bg = "white")
ggsave("figures/5A_new.tiff", p, width = 7.2, height = 4.05, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_5A and figures/5A_new at 600 dpi\n")
