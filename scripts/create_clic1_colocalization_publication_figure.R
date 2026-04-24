#!/usr/bin/env Rscript
# Minimal publication-ready CLIC1 colocalization figure

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})



dir.create("figures", showWarnings = FALSE)
dir.create("SLE_Publication_Package/Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("SLE_Publication_Package/Frontiers_Main_Figures", showWarnings = FALSE, recursive = TRUE)

clean_tissue <- function(x) {
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

validation_path <- "results/targeted_causal_validation.tsv"

if (file.exists(validation_path)) {
  coloc <- fread(validation_path) %>%
    filter(RSID == "rs389884", Gene == "CLIC1") %>%
    transmute(Tissue, PP4)
} else {
  coloc <- fread("results/coloc_results_summary.tsv") %>%
    filter(Locus == "rs389884", Gene == "CLIC1") %>%
    transmute(Tissue, PP4)
}

coloc <- coloc %>%
  mutate(
    Tissue = clean_tissue(Tissue),
    Highlight = PP4 == max(PP4, na.rm = TRUE)
  ) %>%
  arrange(desc(PP4)) %>%
  mutate(Tissue = factor(Tissue, levels = Tissue))

if (nrow(coloc) == 0) {
  stop("No CLIC1 colocalization results found in results/coloc_results_summary.tsv")
}

p <- ggplot(coloc, aes(x = PP4, y = Tissue)) +
  geom_segment(
    aes(x = 0, xend = PP4, y = Tissue, yend = Tissue),
    linewidth = 2.8,
    color = "#d9e6df",
    lineend = "round"
  ) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "#c27d38", linewidth = 0.7) +
  geom_point(
    aes(fill = Highlight),
    shape = 21,
    size = 7.2,
    stroke = 0.6,
    color = "#16324f"
  ) +
  scale_fill_manual(values = c("TRUE" = "#0b7d4b", "FALSE" = "#7fc8a9"), guide = "none") +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "CLIC1 Colocalization",
    x = "Posterior probability of colocalization (PP4)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e8eef3", linewidth = 0.5),
    axis.text.y = element_text(face = "bold", color = "#16324f", size = 12),
    axis.text.x = element_text(color = "#486581"),
    axis.title.x = element_text(color = "#334e68", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 17, color = "#102a43"),
    plot.margin = margin(12, 20, 12, 12)
  )

output_paths <- c(
  "figures/Fig_5_CLIC1_colocalization.png",
  "figures/Fig_5_CLIC1_colocalization.pdf",
  "SLE_Publication_Package/Figures/Fig_5_CLIC1_colocalization.png",
  "SLE_Publication_Package/Figures/Fig_5_CLIC1_colocalization.pdf",
  "SLE_Publication_Package/Frontiers_Main_Figures/Fig_5_CLIC1_colocalization.png",
  "SLE_Publication_Package/Frontiers_Main_Figures/Fig_5_CLIC1_colocalization.pdf"
)

for (path in output_paths) {
  if (grepl("\\.png$", path)) {
    ggsave(path, p, width = 8.2, height = 4.6, dpi = 320, bg = "white")
  } else {
    ggsave(path, p, width = 8.2, height = 4.6, bg = "white")
  }
}

cat("Saved simplified CLIC1 colocalization figure.\n")
