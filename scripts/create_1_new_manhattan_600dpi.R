#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

dir.create("figures", showWarnings = FALSE)

message("Reading full meta-analysis results from zip...")
meta <- fread(
  cmd = "unzip -p GWAS_SLE_summaryStats/Meta_Results.zip Meta_Results.txt",
  select = c("CHR", "SNP", "BP", "P-value"),
  colClasses = list(integer = c("CHR", "BP"), character = "SNP", numeric = "P-value")
)
setnames(meta, c("CHR", "RSID", "BP", "P"))
meta <- meta[!is.na(CHR) & CHR %between% c(1L, 22L) & is.finite(P) & P > 0]

message("Preparing Manhattan coordinates...")
setorder(meta, CHR, BP)
chr_offsets <- meta[, .(chr_len = as.numeric(max(BP, na.rm = TRUE))), by = CHR][order(CHR)]
chr_offsets[, offset := shift(cumsum(chr_len), fill = 0)]
meta <- chr_offsets[, .(CHR, offset)][meta, on = "CHR"]
meta[, bp_cum := as.numeric(BP) + offset]
meta[, logp := -log10(P)]

axis_df <- meta[, .(center = mean(range(bp_cum))), by = CHR]

# Plot all nominally associated variants plus a sparse background, preserving peaks
# while keeping the file light and visually clean.
set.seed(20260525)
background <- meta[P >= 0.05]
if (nrow(background) > 250000) {
  background <- background[sample.int(.N, 250000)]
}
plot_dt <- rbindlist(list(background, meta[P < 0.05]), use.names = TRUE)
plot_dt[, chr_group := factor(CHR %% 2)]

label_source <- fread("results/Supplementary_Table_1.tsv")
label_source <- label_source[is.finite(P_meta)]
label_source[, gene_label := fifelse(
  !is.na(Effector_Gene) & Effector_Gene != "" & Effector_Gene != '""',
  Effector_Gene,
  fifelse(!is.na(Gene) & Gene != "" & Gene != '""', Gene, RSID)
)]
label_source <- label_source[order(P_meta)][1:min(.N, 8)]
label_dt <- merge(
  label_source[, .(RSID, CHR = as.integer(CHR), BP = as.numeric(BP), P = P_meta, gene_label)],
  chr_offsets[, .(CHR, offset)],
  by = "CHR",
  all.x = TRUE
)
label_dt[, `:=`(bp_cum = BP + offset, logp = -log10(P))]

genomewide <- -log10(5e-8)

p <- ggplot(plot_dt, aes(x = bp_cum, y = logp)) +
  geom_point(aes(color = chr_group), alpha = 0.75, size = 0.36, stroke = 0) +
  geom_point(
    data = label_dt,
    aes(x = bp_cum, y = logp),
    inherit.aes = FALSE,
    color = "#111111",
    fill = "#ffffff",
    shape = 21,
    size = 1.05,
    stroke = 0.35
  ) +
  geom_hline(yintercept = genomewide, color = "#d9534f", linetype = "dashed", linewidth = 0.35) +
  geom_text_repel(
    data = label_dt,
    aes(x = bp_cum, y = logp, label = gene_label),
    inherit.aes = FALSE,
    size = 2.45,
    fontface = "bold",
    color = "#111111",
    min.segment.length = 0,
    segment.color = "#6b7280",
    segment.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.22,
    force = 4,
    max.overlaps = Inf,
    seed = 25
  ) +
  scale_color_manual(values = c("0" = "#2f7fb8", "1" = "#2f3f4f"), guide = "none") +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = expansion(mult = c(0.006, 0.006))
  ) +
  scale_y_continuous(
    limits = c(0, 155),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#e8edf3", linewidth = 0.32),
    axis.text.x = element_text(size = 6.8, color = "#334155", face = "bold"),
    axis.text.y = element_text(size = 7.0, color = "#334155", face = "bold"),
    axis.title.x = element_text(size = 8.2, color = "#111111", face = "bold", margin = margin(t = 3)),
    axis.title.y = element_text(size = 8.2, color = "#111111", face = "bold", margin = margin(r = 3)),
    plot.margin = margin(5, 5, 4, 4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("figures/Figure_1.png", p, width = 7.8, height = 4.55, dpi = 600, bg = "white")
ggsave("figures/Figure_1.tiff", p, width = 7.8, height = 4.55, dpi = 600, bg = "white", compression = "lzw")

cat("Saved figures/Figure_1.png and figures/Figure_1.tiff at 600 dpi\n")
