#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(grid)
})

nodes <- data.table(
  name = c(
    "CLIC1", "PTPN22", "IL12A", "HLA-A", "JAK1", "IRF5",
    "STAT1", "TNFR1", "ITGAM", "HLA-B", "TNF", "STAT4",
    "TYK2", "HLA-DR", "IL12B"
  ),
  pathway = c(
    "Cell activation", "JAK-STAT signaling", "Cytokine signaling",
    "Immune response", "JAK-STAT signaling", "Interferon signaling",
    "JAK-STAT signaling", "Cell activation", "Immune response",
    "Immune response", "Cell activation", "JAK-STAT signaling",
    "JAK-STAT signaling", "Immune response", "Cytokine signaling"
  ),
  x = c(0.51, 0.71, 0.83, 0.25, 0.64, 0.83, 0.53, 0.42, 0.36, 0.24, 0.42, 0.49, 0.61, 0.24, 0.83),
  y = c(0.68, 0.22, 0.55, 0.39, 0.38, 0.20, 0.40, 0.17, 0.28, 0.22, 0.33, 0.57, 0.58, 0.58, 0.38)
)

edges <- data.table(
  from = c(
    "CLIC1", "CLIC1", "CLIC1", "CLIC1",
    "PTPN22", "PTPN22", "PTPN22", "PTPN22",
    "IL12A", "IL12A", "IL12A",
    "HLA-A", "HLA-A", "HLA-A",
    "JAK1", "JAK1", "JAK1",
    "IRF5", "IRF5",
    "STAT1", "STAT1", "STAT1", "STAT1", "STAT1", "STAT1", "STAT1", "STAT1",
    "TNFR1",
    "TNF", "TNF", "TNF", "TNF",
    "STAT4", "STAT4",
    "TYK2", "TYK2", "TYK2",
    "HLA-DR",
    "IL12B"
  ),
  to = c(
    "PTPN22", "STAT1", "HLA-A", "JAK1",
    "STAT1", "HLA-A", "IL12A", "JAK1",
    "IL12B", "STAT4", "STAT1",
    "STAT1", "JAK1", "HLA-DR",
    "STAT1", "TYK2", "HLA-A",
    "TNF", "STAT1",
    "STAT4", "TYK2", "TNF", "TNFR1", "ITGAM", "HLA-B", "HLA-A", "IRF5",
    "TNF",
    "STAT1", "STAT4", "ITGAM", "IRF5",
    "TYK2", "IL12B",
    "STAT1", "STAT4", "JAK1",
    "HLA-B",
    "STAT4"
  )
)

graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
V(graph)$degree <- degree(graph)

layout <- nodes[, .(x, y, name)]

pathway_colors <- c(
  "Cell activation" = "#f05b50",
  "Immune response" = "#54c7bd",
  "Cytokine signaling" = "#4db8cf",
  "Interferon signaling" = "#f4bd45",
  "JAK-STAT signaling" = "#6f3fd4"
)

p <- ggraph(graph, layout = "manual", x = layout$x, y = layout$y) +
  geom_edge_link(
    color = "#747b84",
    alpha = 0.48,
    linewidth = 0.56,
    lineend = "round"
  ) +
  geom_node_point(
    aes(fill = pathway, size = degree),
    shape = 21,
    color = "#111111",
    stroke = 0.70,
    alpha = 0.96
  ) +
  geom_node_text(
    aes(label = name),
    size = 2.20,
    fontface = "bold",
    color = "#111111",
    family = "sans"
  ) +
  annotate(
    "text",
    x = 0.51,
    y = 0.745,
    label = "candidate\nsignal",
    size = 2.15,
    lineheight = 0.88,
    fontface = "italic",
    color = "#1f2933",
    family = "sans"
  ) +
  scale_fill_manual(values = pathway_colors, name = NULL) +
  scale_size(range = c(11.2, 15.8), guide = "none") +
  coord_equal(xlim = c(0.06, 0.95), ylim = c(0.075, 0.765), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.key.size = unit(0.14, "in"),
    legend.spacing.y = unit(0.04, "in"),
    legend.text = element_text(size = 7.6, face = "bold", color = "#111111"),
    legend.margin = margin(0, 1, 0, 1),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 4.0)))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_5B.png", p, width = 6.8, height = 4.35, dpi = 600, bg = "white")
ggsave("figures/Figure_5B.tiff", p, width = 6.8, height = 4.35, dpi = 600, bg = "white", compression = "lzw")
ggsave("figures/5B_new.png", p, width = 6.8, height = 4.35, dpi = 600, bg = "white")
ggsave("figures/5B_new.tiff", p, width = 6.8, height = 4.35, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_5B and figures/5B_new at 600 dpi\n")
