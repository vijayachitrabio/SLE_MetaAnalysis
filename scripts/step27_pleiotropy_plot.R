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
})


dir.create("figures", showWarnings = FALSE)

clean_gene <- function(gene, rsid) {
  gene <- as.character(gene)
  gene <- trimws(gsub('^"|"$', "", gene))
  gene[gene %in% c("", "NA", "TBD", "character(0)")] <- NA_character_
  ifelse(is.na(gene), rsid, gene)
}

wrap_trait <- function(x, width = 30) {
  stringr::str_wrap(x, width = width)
}

phewas <- fread("results/phewas_summary_refined.tsv")
master <- fread("results/master_results_table.tsv")

genes_map <- master %>%
  transmute(RSID, Gene = clean_gene(Gene, RSID)) %>%
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
    Trait_Label = wrap_trait(EFO_Trait)
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

p <- ggplot() +
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
  scale_fill_viridis_c(option = "C", end = 0.95, name = expression(-log[10](P))) +
  scale_size_continuous(range = c(3.2, 8.2), guide = "none") +
  labs(
    title = "SLE Genetic Cross-Talk Portfolio",
    subtitle = "GWAS Catalog pleiotropic associations for prioritized SLE loci",
    x = "Prioritized SLE loci",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, face = "bold", color = "#102a43"),
    axis.text.y = element_text(size = 10, color = "#1f2937"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", color = "#102a43", size = 11),
    strip.background = element_rect(fill = "#e5e7eb", color = "#d1d5db"),
    plot.title = element_text(face = "bold", size = 16, color = "#102a43"),
    plot.subtitle = element_text(size = 10.5, color = "#486581"),
    legend.position = c(0.86, 0.16),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  ) +
  guides(fill = guide_colorbar(title.position = "top", barheight = unit(3.2, "cm")))

ggsave("figures/pleiotropy_map_linear.png", p, width = 11, height = 8.2, dpi = 320, bg = "white")
ggsave("figures/pleiotropy_map_linear.pdf", p, width = 11, height = 8.2, bg = "white")

message("Cleaned pleiotropy map saved to figures/pleiotropy_map_linear.png/pdf")
