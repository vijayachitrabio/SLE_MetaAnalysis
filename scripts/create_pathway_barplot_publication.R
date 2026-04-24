#!/usr/bin/env Rscript
# Publication-ready pathway enrichment bar plot

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})


dir.create("figures", showWarnings = FALSE)

wrap_term <- function(x, width = 34) {
  stringr::str_wrap(x, width = width)
}

pathway <- fread("results/pathway_enrichment_filtered.tsv")

plot_df <- pathway %>%
  filter(source %in% c("GO:BP", "REAC")) %>%
  mutate(
    term_plot = term_clean,
    term_plot = str_replace_all(term_plot, "signaling pathway", "signaling"),
    term_plot = str_replace_all(term_plot, "Signalling", "signaling"),
    term_plot = str_replace_all(term_plot, "Immune System", "immune system"),
    term_plot = wrap_term(term_plot),
    source_label = recode(source, "GO:BP" = "GO Biological Process", "REAC" = "Reactome")
  ) %>%
  filter(
    !str_detect(tolower(term_plot), "response to interleukin-12"),
    !str_detect(tolower(term_plot), "cellular response to interleukin-12"),
    !str_detect(tolower(term_plot), "adaptive immune response based on somatic recombinatio")
  ) %>%
  arrange(p_value, desc(intersection_size)) %>%
  group_by(source_label) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  mutate(
    source_label = factor(source_label, levels = c("GO Biological Process", "Reactome"))
  ) %>%
  arrange(log_pval) %>%
  mutate(term_plot = factor(term_plot, levels = term_plot))

source_colors <- c(
  "GO Biological Process" = "#2c7a7b",
  "Reactome" = "#d97706"
)

p <- ggplot(plot_df, aes(x = log_pval, y = term_plot, fill = source_label)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(
    aes(label = intersection_size),
    hjust = -0.18,
    size = 3.4,
    color = "#334155",
    fontface = "bold"
  ) +
  scale_fill_manual(values = source_colors, name = NULL) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(-log[10](P)),
    y = NULL,
    caption = "Numbers at bar ends indicate overlapping SLE-prioritized genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e8edf3", linewidth = 0.5),
    axis.text.y = element_text(color = "#102a43", size = 10),
    axis.text.x = element_text(color = "#475569"),
    axis.title.x = element_text(color = "#334155", margin = margin(t = 8)),
    plot.caption = element_text(size = 9, color = "#64748b"),
    legend.position = "top",
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave("figures/pathway_enrichment_barplot_publication.png", p, width = 10.2, height = 6.8, dpi = 320, bg = "white")
ggsave("figures/pathway_enrichment_barplot_publication.pdf", p, width = 10.2, height = 6.8, bg = "white")

cat("Saved publication-ready pathway bar plot to figures/pathway_enrichment_barplot_publication.png and .pdf\n")
