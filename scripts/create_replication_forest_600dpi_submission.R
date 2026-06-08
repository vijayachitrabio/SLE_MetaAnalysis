#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "SLE_submission"
dir.create(out_dir, showWarnings = FALSE)

dt <- fread("results/Supplementary_Table_5.tsv")
dt <- dt[`Replication status` == "Replicated"]

dt[, Gene := fifelse(grepl("TNFSF4", `Candidate gene/region`), "TNFSF4",
  fifelse(grepl("CLEC16A", `Candidate gene/region`), "CLEC16A",
  fifelse(grepl("IRF8", `Candidate gene/region`), "IRF8",
  fifelse(grepl("IL12A", `Candidate gene/region`), "IL12A",
  fifelse(grepl("ITGAM", `Candidate gene/region`), "ITGAM",
  fifelse(grepl("STAT4", `Candidate gene/region`), "STAT4",
  fifelse(grepl("FCGR2A", `Candidate gene/region`), "FCGR2A", `Candidate gene/region`)))))))]

order_genes <- c("STAT4", "ITGAM", "IRF8", "TNFSF4", "FCGR2A", "IL12A", "CLEC16A")
dt <- dt[Gene %in% order_genes]
dt[, label := paste0(Gene, " (", RSID, ")")]
dt[, label := factor(label, levels = rev(paste0(
  order_genes,
  " (",
  dt[match(order_genes, Gene), RSID],
  ")"
)))]

se_from_beta_p <- function(beta, p) {
  abs(beta / qnorm(p / 2, lower.tail = FALSE))
}

forest_dt <- rbindlist(list(
  dt[, .(
    label,
    Dataset = "Discovery meta-analysis",
    beta = `Beta discovery`,
    se = se_from_beta_p(`Beta discovery`, `P discovery`),
    p_value = `P discovery`
  )],
  dt[, .(
    label,
    Dataset = "Spanish replication",
    beta = `Beta Spanish replication`,
    se = se_from_beta_p(`Beta Spanish replication`, `P Spanish replication`),
    p_value = `P Spanish replication`
  )]
))

forest_dt[, `:=`(
  ci_low = beta - 1.96 * se,
  ci_high = beta + 1.96 * se,
  Strength = fifelse(p_value < 0.05 / 7, "Bonferroni-level", "Nominal")
)]

forest_dt[, Dataset := factor(Dataset, levels = c("Discovery meta-analysis", "Spanish replication"))]
forest_dt[, Strength := factor(Strength, levels = c("Bonferroni-level", "Nominal"))]

p <- ggplot(forest_dt, aes(x = beta, y = label, color = Dataset)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, color = "#7a8492") +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    position = position_dodge(width = 0.55),
    height = 0.14,
    linewidth = 0.72
  ) +
  geom_point(
    aes(shape = Strength),
    position = position_dodge(width = 0.55),
    size = 2.25,
    stroke = 0.75
  ) +
  scale_color_manual(
    values = c(
      "Discovery meta-analysis" = "#1f6fb2",
      "Spanish replication" = "#d16618"
    ),
    name = "Dataset"
  ) +
  scale_shape_manual(
    values = c("Bonferroni-level" = 18, "Nominal" = 16),
    name = "Replication strength"
  ) +
  scale_x_continuous(
    limits = c(-0.58, 0.68),
    breaks = seq(-0.5, 0.5, by = 0.25),
    expand = c(0, 0)
  ) +
  labs(
    title = "Replication forest plot for top reproduced SLE loci",
    x = "Effect size (beta)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 11.2, face = "bold", hjust = 0.5, color = "#111111"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e6eaef", linewidth = 0.35),
    axis.text.y = element_text(size = 7.9, face = "bold", color = "#111827"),
    axis.text.x = element_text(size = 7.4, color = "#4b5563"),
    axis.title.x = element_text(size = 8.2, color = "#111827", face = "bold", margin = margin(t = 6)),
    legend.position = "right",
    legend.title = element_text(size = 7.2, face = "bold", color = "#111827"),
    legend.text = element_text(size = 6.9, color = "#111827"),
    legend.key.size = unit(0.13, "in"),
    legend.spacing.y = unit(0.03, "in"),
    plot.margin = margin(5, 10, 5, 5)
  )

base <- file.path(out_dir, "Figure_4")
ggsave(paste0(base, ".png"), p, width = 6.4, height = 4.25, dpi = 600, bg = "white")
ggsave(paste0(base, ".tiff"), p, width = 6.4, height = 4.25, dpi = 600, bg = "white", compression = "lzw")

cat("Saved 600 dpi replication forest figure to SLE_submission\n")
