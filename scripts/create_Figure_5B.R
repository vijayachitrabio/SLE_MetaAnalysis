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
  x = c(0.60, 0.77, 0.84, 0.77, 0.77, 0.91, 0.58, 0.59, 0.41, 0.24, 0.25, 0.36, 0.43, 0.14, 0.30),
  y = c(0.86, 0.78, 0.72, 0.54, 0.40, 0.31, 0.40, 0.11, 0.16, 0.24, 0.36, 0.44, 0.56, 0.52, 0.79)
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
    color = "#9ca3af",
    alpha = 0.50,
    linewidth = 0.58,
    lineend = "round"
  ) +
  geom_node_point(
    aes(fill = pathway, size = degree),
    shape = 21,
    color = "#111111",
    stroke = 0.75,
    alpha = 0.96
  ) +
  geom_node_text(
    aes(label = name),
    size = 2.35,
    fontface = "bold",
    color = "#111111",
    family = "sans"
  ) +
  scale_fill_manual(values = pathway_colors, name = NULL) +
  scale_size(range = c(13.5, 22), guide = "none") +
  coord_equal(xlim = c(0.10, 0.93), ylim = c(0.09, 0.88), clip = "off") +
  labs(
    title = "Protein-Protein Interaction Network",
    subtitle = "SLE Susceptibility Genes"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5,
      color = "#111111",
      margin = margin(b = 1)
    ),
    plot.subtitle = element_text(
      size = 7.5,
      face = "bold",
      hjust = 0.5,
      color = "#111111",
      margin = margin(b = 7)
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.key.size = unit(0.13, "in"),
    legend.text = element_text(size = 7.5, face = "bold", color = "#111111"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4.2)))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_5B.png", p, width = 6.8, height = 5.8, dpi = 600, bg = "white")
ggsave("figures/Figure_5B.tiff", p, width = 6.8, height = 5.8, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_5B.png and figures/Figure_5B.tiff at 600 dpi\n")
