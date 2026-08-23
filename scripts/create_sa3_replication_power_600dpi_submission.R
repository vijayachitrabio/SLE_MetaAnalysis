#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "SLE_submission"
dir.create(out_dir, showWarnings = FALSE)

dt <- fread("results/Supplementary_Table_5.tsv")
plot_dt <- dt[`Replication status` == "Replicated"]

plot_dt[, Locus := fifelse(
  grepl("TNFSF4", `Candidate gene/region`), "TNFSF4",
  fifelse(grepl("CLEC16A", `Candidate gene/region`), "CLEC16A",
  fifelse(grepl("IRF8", `Candidate gene/region`), "IRF8",
  fifelse(grepl("IL12A", `Candidate gene/region`), "IL12A",
  fifelse(grepl("ITGAM", `Candidate gene/region`), "ITGAM",
  fifelse(grepl("STAT4", `Candidate gene/region`), "STAT4",
  fifelse(grepl("FCGR2A", `Candidate gene/region`), "FCGR2A", `Candidate gene/region`)))))))]

alpha_bonf <- 0.05 / 57
z_bonf <- qnorm(1 - alpha_bonf / 2)
plot_dt[, SE_rep := abs(`Beta Spanish replication` / qnorm(`P Spanish replication` / 2))]
plot_dt[, Statistical_power := (pnorm(abs(`Beta discovery`) / SE_rep - z_bonf) +
  pnorm(-abs(`Beta discovery`) / SE_rep - z_bonf)) * 100]

order_loci <- c("STAT4", "ITGAM", "TNFSF4", "IRF8", "IL12A", "FCGR2A", "CLEC16A")
plot_dt <- plot_dt[Locus %in% order_loci]
plot_dt[, Locus := factor(Locus, levels = rev(order_loci))]
plot_dt[, Power_label := sprintf("%.1f", Statistical_power)]
plot_dt[, label_x := fifelse(Statistical_power > 92, Statistical_power - 3.0, Statistical_power + 2.0)]
plot_dt[, label_hjust := fifelse(Statistical_power > 92, 1, 0)]

p <- ggplot(plot_dt, aes(x = Statistical_power, y = Locus)) +
  geom_col(
    width = 0.70,
    fill = "#4f9656",
    color = NA,
    alpha = 0.98
  ) +
  geom_text(
    aes(x = label_x, label = Power_label, hjust = label_hjust),
    size = 2.75,
    color = "#26323f",
    family = "sans"
  ) +
  scale_x_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    title = "Replication power analysis",
    subtitle = "Power and replication status in the Spanish cohort",
    x = "Statistical power (%)",
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
    legend.position = "none",
    plot.margin = margin(6, 18, 6, 6)
  )

base <- "figures/Supplementary_Figure_4"
ggsave(paste0(base, ".png"), p, width = 4.9, height = 6.2, dpi = 600, bg = "white")
ggsave(paste0(base, ".tiff"), p, width = 4.9, height = 6.2, dpi = 600, bg = "white", compression = "lzw")

cat("Saved SA3 replication power figure to SLE_submission at 600 dpi\n")
