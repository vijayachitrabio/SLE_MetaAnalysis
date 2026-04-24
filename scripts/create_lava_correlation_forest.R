#!/usr/bin/env Rscript
# Forest-style LAVA local genetic correlation figure

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

lava_sum <- fread("results/final_lava_consolidated_loci.tsv")
lava_raw <- fread("results/lava_sle_results.csv")
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

plot_df <- lava_sum %>%
  left_join(master_best, by = "RSID") %>%
  mutate(
    Gene = clean_text(Gene),
    Master_Gene = clean_text(Master_Gene),
    Master_Effector_Gene = clean_text(Master_Effector_Gene),
    Label_Gene = coalesce(Gene, Master_Gene, Master_Effector_Gene, extract_near_gene(Novelty), RSID)
  ) %>%
  left_join(
    lava_raw %>%
      transmute(
        LAVA_LOC = LOC,
        rho = as.numeric(rho),
        rho_lower = as.numeric(`rho.lower`),
        rho_upper = as.numeric(`rho.upper`),
        p = as.numeric(p)
      ),
    by = "LAVA_LOC"
  ) %>%
  mutate(
    P_meta = as.numeric(P_meta),
    LAVA_P = as.numeric(LAVA_P),
    color_value = -log10(pmax(p, 1e-300))
  ) %>%
  filter(!is.na(LAVA_LOC), !is.na(rho), !is.na(rho_lower), !is.na(rho_upper), !is.na(p), rho > 0) %>%
  arrange(LAVA_LOC, P_meta) %>%
  group_by(LAVA_LOC) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(p, desc(rho), P_meta) %>%
  slice_head(n = 10) %>%
  mutate(Label = factor(Label_Gene, levels = rev(Label_Gene)))

if (nrow(plot_df) == 0) {
  stop("No LAVA loci with valid confidence intervals found.")
}

p <- ggplot(plot_df, aes(x = rho, y = Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6b7280", linewidth = 0.7) +
  geom_errorbarh(
    aes(xmin = rho_lower, xmax = rho_upper, color = color_value),
    height = 0,
    linewidth = 1.1
  ) +
  geom_point(aes(color = color_value), size = 3.8) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = expression(-log[10](P))) +
  labs(
    title = "Top LAVA Local Genetic Correlations",
    x = expression(Local~genetic~correlation~(rho)),
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e6ebf0", linewidth = 0.5),
    axis.text.y = element_text(color = "#111827", size = 12),
    axis.text.x = element_text(color = "#4b5563"),
    axis.title.x = element_text(color = "#334155", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 17, color = "#111827"),
    legend.position = c(0.86, 0.18),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(1.6, "cm"),
    plot.margin = margin(12, 16, 12, 12)
  ) +
  guides(color = guide_colorbar(title.position = "top", barwidth = unit(2.4, "cm")))

ggsave("figures/top_lava_correlation_forest.png", p, width = 9.2, height = 6.6, dpi = 320, bg = "white")
ggsave("figures/top_lava_correlation_forest.pdf", p, width = 9.2, height = 6.6, bg = "white")

cat("Saved forest-style LAVA correlation figure to figures/top_lava_correlation_forest.png and .pdf\n")
print(plot_df %>% select(Label_Gene, rho, rho_lower, rho_upper, p))
