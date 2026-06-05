#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dt <- fread("results/Supplementary_Table_6.tsv")
setnames(dt, c("term", "source", "gene_count", "p_value", "fdr"))

wanted <- data.table(
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
  source = c(rep("GO:BP", 6), rep("REAC", 6))
)

plot_dt <- merge(wanted, dt, by = c("term", "source"), sort = FALSE)
plot_dt[, source_label := fifelse(source == "GO:BP", "GO biological process", "Reactome")]
plot_dt[, log_p := -log10(p_value)]
plot_dt[, term_plot := fifelse(
  term == "interleukin-12-mediated signaling pathway",
  "interleukin-12-mediated\nsignaling",
  term
)]
plot_dt[, term_plot := paste0(toupper(substr(term_plot, 1, 1)), substr(term_plot, 2, nchar(term_plot)))]
plot_dt[, term_plot := factor(term_plot, levels = rev(term_plot))]
plot_dt[, source_label := factor(source_label, levels = c("GO biological process", "Reactome"))]

source_colors <- c(
  "GO biological process" = "#2f837a",
  "Reactome" = "#cf741c"
)

p <- ggplot(plot_dt, aes(x = log_p, y = term_plot, fill = source_label)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = gene_count),
    hjust = -0.12,
    size = 3.1,
    color = "#333333"
  ) +
  scale_fill_manual(values = source_colors, name = NULL) +
  scale_x_continuous(
    limits = c(0, 4.25),
    breaks = 0:4,
    expand = c(0, 0)
  ) +
  labs(x = expression(-log[10](P)), y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 9.5) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.key.size = unit(0.14, "in"),
    legend.text = element_text(size = 8, color = "#111111", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e6eaef", linewidth = 0.35),
    axis.text.y = element_text(size = 7.2, color = "#111111", face = "bold"),
    axis.text.x = element_text(size = 7.4, color = "#4b5563", face = "bold"),
    axis.title.x = element_text(size = 8.5, color = "#111111", face = "bold", margin = margin(t = 5)),
    plot.margin = margin(5, 14, 5, 5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_5A.png", p, width = 7.2, height = 4.75, dpi = 600, bg = "white")
ggsave("figures/Figure_5A.tiff", p, width = 7.2, height = 4.75, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_5A.png and figures/Figure_5A.tiff at 600 dpi\n")
