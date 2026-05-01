#!/usr/bin/env Rscript
# scripts/step27_pleiotropy_plot.R
# Cleaner publication-style pleiotropy figure

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(viridis)
  library(cowplot)
})

setwd(".")


dir.create("figures", showWarnings = FALSE)

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub('^"|"$', "", x))
  x[x %in% c("", "NA", "TBD", "character(0)", "Locus_NA")] <- NA_character_
  x
}

extract_near_gene <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- grepl("Known \\(near .+\\)", x)
  out[hit] <- sub("^Known \\(near (.+)\\)$", "\\1", x[hit])
  out
}

choose_gene_label <- function(gene, master_gene, effector_gene, novelty, rsid) {
  label <- dplyr::coalesce(
    clean_text(gene),
    clean_text(master_gene),
    clean_text(effector_gene),
    extract_near_gene(as.character(novelty))
  )
  label <- ifelse(grepl(" / ", label), sub("^.* / ", "", label), label)
  ifelse(is.na(label), rsid, paste0(label, " (", rsid, ")"))
}

wrap_trait <- function(x, width = 30) {
  stringr::str_wrap(x, width = width)
}

capitalize_trait <- function(x) {
  x <- as.character(x)
  idx <- regexpr("[A-Za-z]", x)
  hit <- idx > 0
  if (any(hit)) {
    starts <- substring(x[hit], 1, idx[hit] - 1)
    firsts <- substring(x[hit], idx[hit], idx[hit])
    rests <- substring(x[hit], idx[hit] + 1)
    x[hit] <- paste0(starts, toupper(firsts), rests)
  }
  x
}

phewas <- fread("results/phewas_summary_refined.tsv")
master <- fread("results/master_results_table.tsv")

genes_map <- master %>%
  mutate(
    Master_Gene = clean_text(Gene),
    Master_Effector_Gene = clean_text(Effector_Gene),
    has_effector = !is.na(Master_Effector_Gene),
    has_gene = !is.na(Master_Gene)
  ) %>%
  arrange(desc(has_effector), desc(has_gene)) %>%
  group_by(RSID) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(RSID, Gene = choose_gene_label(Gene, Master_Gene, Master_Effector_Gene, Novelty, RSID)) %>%
  distinct()

whitelist_keywords <- c(
  "lupus", "arthritis", "sclerosis", "sjogren", "vitiligo", "psoriasis",
  "crohn", "colitis", "thyroid", "kawasaki", "diabetes",
  "leukocyte", "neutrophil", "lymphocyte", "monocyte", "platelet",
  "hematology", "glomerulonephritis", "sarcoidosis", "inflammatory bowel"
)

blacklist_keywords <- c(
  "protein amount", "protein level", "ratio", "lipids", "fatty acids",
  "inguinal hernia", "medication", "takes medication"
)

plot_data <- phewas %>%
  inner_join(genes_map, by = "RSID") %>%
  filter(grepl(paste(whitelist_keywords, collapse = "|"), EFO_Trait, ignore.case = TRUE)) %>%
  filter(!grepl(paste(blacklist_keywords, collapse = "|"), EFO_Trait, ignore.case = TRUE)) %>%
  filter(is.finite(P_value), P_value < 1e-6) %>%
  mutate(
    logP = pmin(-log10(P_value), 60),
    Category = ifelse(Category == "Immune-Mediated", "Immune-mediated", "Other trait"),
    Trait_Label = wrap_trait(capitalize_trait(EFO_Trait))
  ) %>%
  group_by(Gene, EFO_Trait, Trait_Label, Category) %>%
  summarise(logP = max(logP, na.rm = TRUE), .groups = "drop")

top_traits <- plot_data %>%
  group_by(EFO_Trait, Trait_Label, Category) %>%
  summarise(
    peak_logP = max(logP, na.rm = TRUE),
    n_loci = n_distinct(Gene),
    .groups = "drop"
  ) %>%
  arrange(desc(n_loci), desc(peak_logP)) %>%
  slice_head(n = 12)

plot_data_sub <- plot_data %>%
  semi_join(top_traits, by = c("EFO_Trait", "Trait_Label", "Category"))

gene_order <- plot_data_sub %>%
  group_by(Gene) %>%
  summarise(peak_logP = max(logP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(peak_logP)) %>%
  pull(Gene)

trait_order <- top_traits %>%
  arrange(Category, peak_logP) %>%
  pull(Trait_Label)

plot_data_sub <- plot_data_sub %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Trait_Label = factor(Trait_Label, levels = trait_order)
  )

lane_df <- expand.grid(
  Gene = levels(plot_data_sub$Gene),
  Trait_Label = levels(plot_data_sub$Trait_Label),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(top_traits %>% select(Trait_Label, Category), by = "Trait_Label") %>%
  mutate(
    Gene = factor(Gene, levels = levels(plot_data_sub$Gene)),
    Trait_Label = factor(Trait_Label, levels = levels(plot_data_sub$Trait_Label))
  )

p_base <- ggplot() +
  geom_point(
    data = lane_df,
    aes(x = Gene, y = Trait_Label),
    shape = 21,
    size = 1.7,
    stroke = 0.15,
    fill = "#eef2f7",
    color = "#d9e2ec"
  ) +
  geom_point(
    data = plot_data_sub,
    aes(x = Gene, y = Trait_Label, size = logP, fill = logP),
    shape = 21,
    color = "#1f2937",
    stroke = 0.35,
    alpha = 0.98
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_viridis_c(option = "C", end = 0.95, name = NULL) +
  scale_size_continuous(range = c(3.2, 8.2), guide = "none") +
  labs(
    title = NULL,
    x = "Prioritized SLE loci",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 34,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "#102a43",
      margin = margin(t = 6)
    ),
    axis.title.x = element_text(color = "#102a43", margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, color = "#1f2937"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", color = "#102a43", size = 11),
    strip.background = element_rect(fill = "#e5e7eb", color = "#d1d5db"),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    legend.margin = margin(b = 6),
    plot.margin = margin(6, 14, 16, 10)
  ) +
  guides(fill = guide_colorbar(barwidth = unit(4.8, "cm"), barheight = unit(0.35, "cm")))

legend_grob <- cowplot::get_legend(
  p_base +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      plot.margin = margin(0, 0, 0, 0)
    )
)

legend_row <- cowplot::ggdraw() +
  cowplot::draw_grob(legend_grob, x = 0, y = 1, width = 0.42, height = 1, hjust = 0, vjust = 1)

p_main <- p_base + theme(legend.position = "none")

p_final <- cowplot::plot_grid(
  legend_row,
  p_main,
  ncol = 1,
  align = "v",
  rel_heights = c(0.08, 1)
)

ggsave("figures/pleiotropy_map_linear.png", p_final, width = 11, height = 8.8, dpi = 320, bg = "white")
ggsave("figures/pleiotropy_map_linear.pdf", p_final, width = 11, height = 8.8, bg = "white")

message("Cleaned pleiotropy map saved to figures/pleiotropy_map_linear.png/pdf")
