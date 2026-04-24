#!/usr/bin/env Rscript
# Clean top LAVA-supported loci figure

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})



dir.create("figures", showWarnings = FALSE)

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub('^"|"$', "", x))
  x[x %in% c("", "NA", "TBD", "character(0)", "Locus_NA")] <- NA_character_
  x
}

extract_near_gene <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  hit <- grepl("Known \\(near .+\\)", x)
  out[hit] <- sub("^Known \\(near (.+)\\)$", "\\1", x[hit])
  out
}

lava <- fread("results/final_lava_consolidated_loci.tsv")
master <- fread("results/master_results_table.tsv")

master_best <- master %>%
  mutate(
    Gene = clean_text(Gene),
    Effector_Gene = clean_text(Effector_Gene),
    has_effector = !is.na(Effector_Gene),
    has_gene = !is.na(Gene)
  ) %>%
  arrange(desc(has_effector), desc(has_gene)) %>%
  group_by(RSID) %>%
  slice(1) %>%
  ungroup() %>%
  select(RSID, Master_Gene = Gene, Master_Effector_Gene = Effector_Gene)

plot_df <- lava %>%
  left_join(master_best, by = "RSID") %>%
  mutate(
    Gene = clean_text(Gene),
    Master_Gene = clean_text(Master_Gene),
    Master_Effector_Gene = clean_text(Master_Effector_Gene),
    Label_Gene = coalesce(Gene, Master_Gene, Master_Effector_Gene, extract_near_gene(Novelty), RSID),
    Label = paste0(Label_Gene, "  ", RSID),
    LAVA_P = as.numeric(LAVA_P),
    LAVA_rg = as.numeric(LAVA_rg),
    P_meta = as.numeric(P_meta),
    Novelty_Group = ifelse(Novelty == "Putative Novel", "Putative novel", "Known / nearby known")
  ) %>%
  filter(!is.na(LAVA_LOC), !is.na(LAVA_P), !is.na(LAVA_rg), LAVA_rg > 0) %>%
  arrange(LAVA_LOC, P_meta) %>%
  group_by(LAVA_LOC) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(LAVA_P, desc(LAVA_rg), P_meta) %>%
  slice_head(n = 10) %>%
  mutate(
    minus_log10_lava_p = -log10(LAVA_P),
    Label = factor(Label, levels = rev(Label))
  )

if (nrow(plot_df) == 0) {
  stop("No positive LAVA loci with valid P-values found.")
}

bonf_threshold <- -log10(0.05 / nrow(lava))

p <- ggplot(plot_df, aes(x = minus_log10_lava_p, y = Label)) +
  geom_segment(
    aes(x = 0, xend = minus_log10_lava_p, y = Label, yend = Label),
    linewidth = 0.9,
    color = "#d8dee6"
  ) +
  geom_vline(xintercept = bonf_threshold, linetype = "dashed", color = "#c27d38", linewidth = 0.7) +
  geom_point(
    aes(fill = Novelty_Group),
    shape = 21,
    size = 6.2,
    stroke = 0.55,
    color = "#1f2d3d"
  ) +
  scale_fill_manual(values = c(
    "Putative novel" = "#d95f02",
    "Known / nearby known" = "#1b9e77"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "Top LAVA-Supported Loci",
    x = expression(-log[10](LAVA~P)),
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#edf1f5", linewidth = 0.5),
    axis.text.y = element_text(face = "bold", color = "#18212b"),
    axis.text.x = element_text(color = "#4b5563"),
    axis.title.x = element_text(color = "#334155", margin = margin(t = 8)),
    plot.title = element_text(face = "bold", size = 16, color = "#102a43"),
    legend.position = "top",
    plot.margin = margin(12, 16, 12, 12)
  )

ggsave("figures/top_lava_loci.png", p, width = 9.2, height = 6.4, dpi = 320, bg = "white")
ggsave("figures/top_lava_loci.pdf", p, width = 9.2, height = 6.4, bg = "white")

cat("Saved top LAVA loci figure to figures/top_lava_loci.png and .pdf\n")
print(plot_df %>% select(RSID, Label_Gene, LAVA_LOC, LAVA_rg, LAVA_P, Novelty))
