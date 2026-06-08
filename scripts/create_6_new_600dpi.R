#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

phewas <- fread("results/Supplementary_Table_7.tsv")

loci <- data.table(
  RSID = c(
    "rs389884", "rs34572943", "rs17849501", "rs4853458",
    "rs10912578", "rs13332649", "rs41272536", "rs6671847"
  ),
  locus_label = c(
    "CLIC1 (rs389884)",
    "ITGAM (rs34572943)",
    "rs17849501",
    "STAT4 (rs4853458)",
    "TNFSF4 (rs10912578)",
    "IRF8 (rs13332649)",
    "rs41272536",
    "FCGR2A (rs6671847)"
  )
)

traits <- data.table(
  EFO_Trait = c(
    "systemic lupus erythematosus",
    "cutaneous lupus erythematosus",
    "rheumatoid arthritis",
    "ACPA-positive rheumatoid arthritis",
    "membranous glomerulonephritis",
    "systemic sclerosis",
    "sarcoidosis",
    "hypothyroidism",
    "thyroid gland disorder",
    "neutrophil count",
    "myeloid leukocyte count",
    "leukocyte quantity"
  ),
  trait_label = c(
    "Systemic lupus erythematosus",
    "Cutaneous lupus erythematosus",
    "Rheumatoid arthritis",
    "ACPA-positive rheumatoid\narthritis",
    "Membranous glomerulonephritis",
    "Systemic sclerosis",
    "Sarcoidosis",
    "Hypothyroidism",
    "Thyroid gland disorder",
    "Neutrophil count",
    "Myeloid leukocyte count",
    "Leukocyte quantity"
  ),
  category = c(rep("Immune-mediated", 4), rep("Other trait", 8))
)

plot_dt <- merge(phewas, loci, by = "RSID")
plot_dt <- merge(plot_dt, traits, by = "EFO_Trait")
plot_dt[, log_p := pmin(-log10(P_value), 60)]
plot_dt <- plot_dt[, .(log_p = max(log_p, na.rm = TRUE)), by = .(RSID, locus_label, EFO_Trait, trait_label, category)]

grid_dt <- CJ(locus_label = loci$locus_label, trait_label = traits$trait_label, unique = TRUE)
grid_dt <- merge(grid_dt, loci[, .(locus_label)], by = "locus_label", all.x = TRUE)
grid_dt <- merge(grid_dt, traits[, .(trait_label, category)], by = "trait_label", all.x = TRUE)
grid_dt[, x := match(locus_label, loci$locus_label)]
grid_dt[, y := match(trait_label, rev(traits$trait_label))]

plot_dt[, x := match(locus_label, loci$locus_label)]
plot_dt[, y := match(trait_label, rev(traits$trait_label))]

cat_dt <- traits[, .(
  ymin = min(match(trait_label, rev(traits$trait_label))) - 0.45,
  ymax = max(match(trait_label, rev(traits$trait_label))) + 0.45
), by = category]
cat_dt[, x_label := -1.32]
cat_dt[, x_min := -2.18]
cat_dt[, x_max := -0.34]

p <- ggplot() +
  geom_rect(
    data = cat_dt,
    aes(xmin = x_min, xmax = x_max, ymin = ymin, ymax = ymax),
    fill = "#e5e7eb",
    color = "#d1d5db",
    linewidth = 0.35
  ) +
  geom_text(
    data = cat_dt,
    aes(x = x_label, y = (ymin + ymax) / 2, label = category),
    color = "#102a43",
    fontface = "bold",
    size = 2.35
  ) +
  geom_point(
    data = grid_dt,
    aes(x = x, y = y),
    shape = 21,
    size = 2.1,
    stroke = 0.20,
    fill = "#eef3f8",
    color = "#dbe5ef"
  ) +
  geom_point(
    data = plot_dt,
    aes(x = x, y = y, size = log_p, fill = log_p),
    shape = 21,
    color = "#222222",
    stroke = 0.45,
    alpha = 0.98
  ) +
  scale_fill_viridis_c(
    option = "C",
    end = 0.96,
    limits = c(6, 60),
    breaks = c(10, 20, 30, 40, 50, 60),
    name = expression(-log[10](P))
  ) +
  scale_size_continuous(
    limits = c(6, 60),
    range = c(2.8, 7.8),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq_len(nrow(loci)),
    labels = loci$locus_label,
    limits = c(-2.25, nrow(loci) + 0.45),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(nrow(traits)),
    labels = rev(traits$trait_label),
    limits = c(0.45, nrow(traits) + 0.55),
    expand = c(0, 0)
  ) +
  labs(x = "Prioritized SLE loci", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      color = "#102a43",
      size = 6.8,
      angle = 42,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(color = "#1f2937", size = 7.1),
    axis.title.x = element_text(face = "bold", color = "#102a43", size = 8.2, margin = margin(t = 3)),
    axis.ticks = element_blank(),
    legend.position = c(0.01, 1.105),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
    legend.title = element_text(size = 6.8, color = "#111111", face = "bold"),
    legend.text = element_text(size = 6.2, color = "#111111", face = "bold"),
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(19, 7, 4, 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    barwidth = unit(1.45, "in"),
    barheight = unit(0.10, "in"),
    ticks.colour = "#ffffff",
    frame.colour = NA
  ))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_6.png", p, width = 7.25, height = 5.35, dpi = 600, bg = "white")
ggsave("figures/Figure_6.tiff", p, width = 7.25, height = 5.35, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_6.png and figures/Figure_6.tiff at 600 dpi\n")
