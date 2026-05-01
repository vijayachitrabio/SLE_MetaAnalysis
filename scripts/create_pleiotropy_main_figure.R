#!/usr/bin/env Rscript
# Compact main-text pleiotropy figure focused on the strongest immune-mediated overlaps.

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

# Main-text version: keep only the clearest immune-mediated overlaps and named loci.
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
  )

p <- ggplot() +
  geom_point(
    data = grid_df,
    aes(x = Gene, y = Trait_Label),
    shape = 21,
    size = 2.4,
    stroke = 0.18,
    fill = "#edf2f7",
    color = "#d9e2ec"
  ) +
  geom_point(
    data = plot_data_focus,
    aes(x = Gene, y = Trait_Label, size = logP, fill = logP),
    shape = 21,
    color = "#1f2937",
    stroke = 0.5,
    alpha = 0.98
  ) +
  scale_fill_viridis_c(
    option = "C",
    end = 0.96,
    limits = c(6, 60),
    breaks = c(10, 20, 30, 40, 50, 60),
    name = expression(-log[10](P))
  ) +
  scale_size_continuous(
    limits = c(6, 60),
    range = c(4.3, 13.5),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", color = "#102a43", size = 11),
    axis.text.y = element_text(color = "#1f2937", size = 10.5),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, color = "#334155"),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(6, 12, 8, 6)
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(5.2, "cm"), barheight = unit(0.35, "cm")))

ggsave("figures/pleiotropy_main_version.png", p, width = 7.6, height = 4.9, dpi = 320, bg = "white")
ggsave("figures/pleiotropy_main_version.pdf", p, width = 7.6, height = 4.9, bg = "white")

message("Saved compact main-text pleiotropy figure to figures/pleiotropy_main_version.png/pdf")
