#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

base_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis"
out_dir <- file.path(base_dir, "2026-05-18_figure_edits")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sentence_case_term <- function(x) {
  x <- str_replace_all(x, "[Ss]ignalling", "signaling")
  x <- str_replace_all(x, "signaling pathway", "signaling")
  x <- str_replace_all(x, "Immune System", "immune system")
  x <- tolower(x)
  x <- str_replace_all(x, "\\bgo\\b", "GO")
  x <- str_replace_all(x, "\\bsle\\b", "SLE")
  x <- str_replace_all(x, "\\bjak-stat\\b", "JAK-STAT")
  x
}

pathway <- fread(file.path(base_dir, "results", "pathway_enrichment_filtered.tsv"))

plot_df <- pathway %>%
  filter(source %in% c("GO:BP", "REAC")) %>%
  mutate(
    term_plot = sentence_case_term(term_clean),
    source_label = recode(source, "GO:BP" = "GO biological process", "REAC" = "Reactome")
  ) %>%
  filter(
    !str_detect(tolower(term_plot), "response to interleukin-12$"),
    !str_detect(tolower(term_plot), "cellular response to interleukin-12$"),
    !str_detect(tolower(term_plot), "adaptive immune response based on somatic")
  ) %>%
  arrange(p_value, desc(intersection_size)) %>%
  group_by(source_label) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  mutate(
    source_label = factor(source_label, levels = c("GO biological process", "Reactome")),
    term_plot = str_wrap(term_plot, width = 31)
  ) %>%
  arrange(source_label, desc(log_pval)) %>%
  mutate(term_plot = factor(term_plot, levels = rev(term_plot)))

source_colors <- c(
  "GO biological process" = "#2f7f7a",
  "Reactome" = "#c86f1d"
)

plot_pathway_bar <- ggplot(plot_df, aes(x = log_pval, y = term_plot, fill = source_label)) +
  geom_col(width = 0.68, color = NA, alpha = 0.96) +
  geom_text(
    aes(label = intersection_size),
    hjust = -0.18,
    size = 2.85,
    color = "#344054"
  ) +
  scale_fill_manual(values = source_colors, name = NULL) +
  scale_x_continuous(
    expression(-log[10](italic(P))),
    limits = c(0, max(plot_df$log_pval) + 0.28),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e6eaf0", linewidth = 0.35),
    axis.text.y = element_text(size = 7.8, color = "#263238", lineheight = 0.95),
    axis.text.x = element_text(size = 8, color = "#4b5563"),
    axis.title.x = element_text(size = 8.6, color = "#374151", margin = margin(t = 8)),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 7.6, color = "#263238"),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  )

png_path <- file.path(out_dir, "Pathway_enrichment_bar_minor_improvements_20260518.png")
pdf_path <- file.path(out_dir, "Pathway_enrichment_bar_minor_improvements_20260518.pdf")

ggsave(png_path, plot_pathway_bar, width = 7.2, height = 4.8, dpi = 450, bg = "white")
ggsave(pdf_path, plot_pathway_bar, width = 7.2, height = 4.8, bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
