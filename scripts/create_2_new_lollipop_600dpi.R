#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub('^"|"$', "", x))
  x[x %in% c("", "NA", "TBD", "character(0)", "Locus_NA")] <- NA_character_
  x
}

dt <- fread("results/Supplementary_Table_1.tsv")
dt[, `:=`(
  P_meta = as.numeric(P_meta),
  BETA = as.numeric(BETA),
  Gene = clean_text(Gene),
  Effector_Gene = clean_text(Effector_Gene),
  Novelty = clean_text(Novelty),
  LAVA_Status = clean_text(LAVA_Status)
)]
dt[, plot_gene := fifelse(!is.na(Effector_Gene), Effector_Gene, fifelse(!is.na(Gene), Gene, RSID))]
dt[, label := paste(plot_gene, RSID)]
dt[, log_p := -log10(P_meta)]
dt[, abs_beta := abs(BETA)]
dt[, locus_class := fifelse(Novelty == "Putative Novel", "Putative novel", "Known / nearby known")]
dt[, support := fifelse(grepl("Confirmed", LAVA_Status), "LAVA replicated", "Not replicated")]

plot_dt <- dt[order(P_meta)][1:10]
plot_dt[, label := factor(label, levels = rev(label))]
plot_dt[, beta_label_x := log_p + 4.2 + (abs_beta * 4.2)]
x_limit <- max(plot_dt$beta_label_x) + 7

class_colors <- c(
  "Known / nearby known" = "#1f9d7a",
  "Putative novel" = "#d26f1d"
)
support_shapes <- c(
  "LAVA replicated" = 21,
  "Not replicated" = 24
)

p <- ggplot(plot_dt, aes(x = log_p, y = label)) +
  geom_segment(
    aes(x = 0, xend = log_p, y = label, yend = label, color = locus_class),
    linewidth = 0.72,
    alpha = 0.88
  ) +
  geom_point(
    aes(fill = locus_class, shape = support, size = abs_beta),
    color = "#111111",
    stroke = 0.45,
    alpha = 0.98
  ) +
  geom_text(
    aes(x = beta_label_x, label = sprintf("%.2f", abs_beta)),
    size = 2.35,
    color = "#333333",
    fontface = "bold",
    hjust = 0
  ) +
  scale_color_manual(values = class_colors, guide = "none") +
  scale_fill_manual(values = class_colors, name = "Locus class") +
  scale_shape_manual(values = support_shapes, name = "Support") +
  scale_size_continuous(name = expression(paste("|", beta, "|")), range = c(3.4, 6.7), breaks = c(0.4, 0.6, 0.8)) +
  scale_x_continuous(
    breaks = seq(0, 160, by = 40),
    limits = c(0, x_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Top 10 SLE loci ranked by statistical strength",
    x = expression(-log[10](P[meta])),
    y = NULL
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e8edf3", linewidth = 0.32),
    axis.text.y = element_text(size = 6.4, color = "#111111", face = "bold"),
    axis.text.x = element_text(size = 6.5, color = "#334155", face = "bold"),
    axis.title.x = element_text(size = 7.4, color = "#111111", face = "bold", margin = margin(t = 4)),
    plot.title = element_text(size = 8.5, color = "#111111", face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 6.4, color = "#111111", face = "bold"),
    legend.text = element_text(size = 5.8, color = "#111111"),
    legend.key.size = unit(0.15, "in"),
    legend.key.height = unit(0.17, "in"),
    legend.spacing.y = unit(0.04, "in"),
    legend.box.spacing = unit(0.10, "in"),
    plot.margin = margin(5, 5, 4, 4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 21, size = 3.2, color = "#111111", stroke = 0.45)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(size = 3.2, fill = "#7b8491", color = "#111111", stroke = 0.45)
    ),
    size = guide_legend(
      order = 3,
      override.aes = list(shape = 21, fill = "#111111", color = "#111111", stroke = 0.25)
    )
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_2.png", p, width = 7.25, height = 4.6, dpi = 600, bg = "white")
ggsave("figures/Figure_2.tiff", p, width = 7.25, height = 4.6, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_2.png and figures/Figure_2.tiff at 600 dpi\n")
