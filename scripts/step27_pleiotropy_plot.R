#!/usr/bin/env Rscript
# scripts/step27_pleiotropy_plot.R
# Visualizing pleiotropic associations of SLE loci with other immune/blood traits

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# 1. Load Data
phewas <- fread("results/phewas_summary_refined.tsv")
master <- fread("results/master_results_table.tsv")

# Identify High-Confidence Genes
genes_map <- master %>%
  select(RSID, Gene) %>%
  mutate(Gene = ifelse(Gene == "" | is.na(Gene), RSID, Gene)) %>%
  unique()

# 2. Filter & Refine Traits
# Only include Immune-Mediated, Blood markers, or relevant inflammatory signals
whitelist_keywords <- c("Lupus", "Arthritis", "Sclerosis", "Thyroid", "Diabetes", "Crohn", 
                        "Colitis", "Psoriasis", "Leukocyte", "Neutrophil", "Lymphocyte", 
                        "Monocyte", "Platelet", "Hematology", "C-reactive", "Chemokine", 
                        "Interferon", "Autoimmune", "Ankylosing", "Sjogren", "Vitiligo",
                        "protein level", "quantity", "count")

plot_data <- phewas %>%
  inner_join(genes_map, by = "RSID") %>%
  filter(grepl(paste(whitelist_keywords, collapse="|"), EFO_Trait, ignore.case = TRUE)) %>%
  # Remove broad/vague terms
  filter(!grepl("hernia|hearing|body mass|stature|intelligence", EFO_Trait, ignore.case=TRUE)) %>%
  # Filter for significance
  filter(P_value < 1e-6) %>%
  mutate(logP = -log10(P_value)) %>%
  mutate(logP = ifelse(logP > 100, 100, logP)) # Cap for visual scale

# Top traits only to avoid clutter (Top 20 most frequent traits)
top_traits <- plot_data %>%
  group_by(EFO_Trait) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  head(20) %>%
  pull(EFO_Trait)

plot_data_sub <- plot_data %>%
  filter(EFO_Trait %in% top_traits)

# 3. Visualization
message("Generating Pleiotropy Bubble Plot...")

p <- ggplot(plot_data_sub, aes(x = Gene, y = EFO_Trait)) +
  geom_point(aes(size = logP, color = Category), alpha = 0.7) +
  scale_size_continuous(range = c(2, 10), breaks = c(10, 30, 50, 80)) +
  scale_color_manual(values = c("Immune-Mediated" = "#e74c3c", "Other Trait" = "#3498db")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    title = "Trait Correlation Across High-Confidence SLE Loci",
    subtitle = "Pleiotropic associations from GWAS Catalog (P < 1e-6)",
    x = "SLE Risk Loci (Gene Symbol)",
    y = "Associated Trait (EFO)",
    size = "-log10(P-value)",
    color = "Trait Category"
  )

# Save
dir.create("figures", showWarnings = FALSE)
ggsave("figures/pleiotropy_map.png", p, width = 12, height = 9, dpi = 300)
ggsave("figures/pleiotropy_map.pdf", p, width = 12, height = 9)

message("Pleiotropy Map saved to figures/pleiotropy_map.png")
