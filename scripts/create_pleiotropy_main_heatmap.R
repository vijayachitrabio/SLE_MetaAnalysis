#!/usr/bin/env Rscript
# Compact main-text pleiotropy figure as a heatmap rather than a dot plot.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
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
  ifelse(is.na(label), rsid, label)
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
    Category = ifelse(Category == "Immune-Mediated", "Immune-mediated", "Other trait")
  ) %>%
  group_by(Gene, RSID, EFO_Trait, Category) %>%
  summarise(logP = max(logP, na.rm = TRUE), .groups = "drop")

focus_traits <- c(
  "systemic lupus erythematosus",
  "cutaneous lupus erythematosus",
  "rheumatoid arthritis",
  "ACPA-positive rheumatoid arthritis",
  "ulcerative colitis"
)

focus_genes <- c("ITGAM", "STAT4", "TNFSF4", "IRF8", "FCGR2A")

plot_data_focus <- plot_data %>%
  filter(Category == "Immune-mediated") %>%
  filter(EFO_Trait %in% focus_traits, Gene %in% focus_genes) %>%
  mutate(
    Trait_Label = recode(
      EFO_Trait,
      "systemic lupus erythematosus" = "Systemic lupus erythematosus",
      "cutaneous lupus erythematosus" = "Cutaneous lupus erythematosus",
      "rheumatoid arthritis" = "Rheumatoid arthritis",
      "ACPA-positive rheumatoid arthritis" = "ACPA-positive rheumatoid arthritis",
      "ulcerative colitis" = "Ulcerative colitis"
    ),
    Gene = factor(Gene, levels = focus_genes),
    Trait_Label = factor(Trait_Label, levels = rev(c(
      "Systemic lupus erythematosus",
      "Cutaneous lupus erythematosus",
      "Rheumatoid arthritis",
      "ACPA-positive rheumatoid arthritis",
      "Ulcerative colitis"
    )))
  )

grid_df <- expand.grid(
  Gene = focus_genes,
  Trait_Label = levels(plot_data_focus$Trait_Label),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    Gene = factor(Gene, levels = focus_genes),
    Trait_Label = factor(Trait_Label, levels = levels(plot_data_focus$Trait_Label))
  ) %>%
  left_join(plot_data_focus %>% select(Gene, Trait_Label, logP), by = c("Gene", "Trait_Label")) %>%
  mutate(
    present = !is.na(logP),
    label = ifelse(present, sprintf("%.1f", logP), "")
  )

p <- ggplot(grid_df, aes(x = Gene, y = Trait_Label)) +
  geom_tile(fill = "#f3f6fa", color = "#d9e2ec", linewidth = 0.9, width = 0.92, height = 0.92) +
  geom_tile(
    data = subset(grid_df, present),
    aes(fill = logP),
    color = "#243447",
    linewidth = 1.0,
    width = 0.92,
    height = 0.92
  ) +
  geom_text(
    data = subset(grid_df, present),
    aes(label = label, color = logP > 34),
    size = 4.8,
    fontface = "bold"
  ) +
  scale_fill_viridis_c(
    option = "E",
    begin = 0.08,
    end = 0.92,
    limits = c(6, 60),
    breaks = c(10, 20, 30, 40, 50, 60),
    name = expression(-log[10](P))
  ) +
  scale_color_manual(values = c("FALSE" = "#ffffff", "TRUE" = "#102a43"), guide = "none") +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "#102a43", size = 11),
    axis.text.y = element_text(color = "#1f2937", size = 10.5),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification = "right",
    legend.box.just = "right",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, color = "#334155"),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    legend.margin = margin(t = 4),
    plot.margin = margin(0, 12, 2, 6)
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(4.4, "cm"), barheight = unit(0.28, "cm")))

ggsave("figures/pleiotropy_main_heatmap.png", p, width = 7.2, height = 4.4, dpi = 320, bg = "white")
ggsave("figures/pleiotropy_main_heatmap.pdf", p, width = 7.2, height = 4.4, bg = "white")

message("Saved main-text pleiotropy heatmap to figures/pleiotropy_main_heatmap.png/pdf")
