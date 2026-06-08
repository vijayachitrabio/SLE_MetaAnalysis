#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "SLE_submission"
dir.create(out_dir, showWarnings = FALSE)

dt <- fread("results/Supplementary_Table_2.tsv")
setnames(dt, "I²", "I2")

gene_map <- data.table(
  RSID = c(
    "rs4853458", "rs34572943", "rs13332649", "rs10912578",
    "rs6671847", "rs2647928", "rs12928726"
  ),
  Locus = c("STAT4", "ITGAM", "IRF8", "TNFSF4", "FCGR2A", "IL12A", "CLEC16A"),
  Display_order = 1:7
)

plot_dt <- merge(gene_map, dt, by = "RSID", all.x = TRUE, sort = FALSE)
setorder(plot_dt, Display_order)
plot_dt[, Locus := factor(Locus, levels = rev(Locus))]
plot_dt[, I2_label := ifelse(I2 > 0, sprintf("%.1f", I2), "0.0")]
plot_dt[, label_x := fifelse(I2 >= 85, I2 + 1.8, fifelse(I2 > 60, I2 - 2.2, 1.8))]
plot_dt[, label_hjust := fifelse(I2 > 60 & I2 < 85, 1, 0)]

p <- ggplot(plot_dt, aes(x = I2, y = Locus)) +
  geom_col(width = 0.70, fill = "#6f9fc6", color = NA, alpha = 0.98) +
  geom_vline(
    xintercept = c(50, 75),
    linetype = "dashed",
    linewidth = 0.55,
    color = "#cc2f32"
  ) +
  annotate(
    "text",
    x = 50,
    y = 7.28,
    label = "Moderate (50%)",
    size = 2.55,
    color = "#8b1e24",
    family = "sans"
  ) +
  annotate(
    "text",
    x = 75,
    y = 7.28,
    label = "High (75%)",
    size = 2.55,
    color = "#8b1e24",
    family = "sans"
  ) +
  geom_text(
    aes(x = label_x, label = I2_label, hjust = label_hjust),
    size = 2.75,
    color = "#26323f",
    family = "sans"
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    title = "SA1: Heterogeneity by locus",
    subtitle = "Fixed-effects discovery meta-analysis",
    x = expression(I^2~"(%)"),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5,
      color = "#111111",
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 8.5,
      hjust = 0.02,
      color = "#2f3a45",
      margin = margin(b = 8)
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e7ebef", linewidth = 0.35),
    axis.text.y = element_text(size = 8.3, color = "#1f2933", face = "bold"),
    axis.text.x = element_text(size = 7.6, color = "#4b5563"),
    axis.title.x = element_text(size = 8.6, color = "#111111", face = "bold", margin = margin(t = 6)),
    plot.margin = margin(6, 18, 6, 6)
  )

base <- "figures/Supplementary_Figure_3"
ggsave(paste0(base, ".png"), p, width = 6.2, height = 4.2, dpi = 600, bg = "white")
ggsave(paste0(base, ".tiff"), p, width = 6.2, height = 4.2, dpi = 600, bg = "white", compression = "lzw")
cat("Saved SA1 heterogeneity figure to SLE_submission at 600 dpi\n")
