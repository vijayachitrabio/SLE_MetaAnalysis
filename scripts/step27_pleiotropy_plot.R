#!/usr/bin/env Rscript
# scripts/step27_pleiotropy_plot.R
# Final Linear Pleiotropy Map (High-aesthetic Bubble Plot with Horizontal Traits)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
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
whitelist_keywords <- c("Lupus", "SLE", "Arthritis", "Sclerosis", "Thyroid", "Diabetes", "Crohn", 
                        "Colitis", "Psoriasis", "Leukocyte", "Neutrophil", "Lymphocyte", 
                        "Monocyte", "Platelet", "Hematology", "C-reactive", "Chemokine", 
                        "Interferon", "Autoimmune", "Ankylosing", "Sjogren", "Vitiligo",
                        "membranous glomerulonephritis", "Scleroderma")

plot_data <- phewas %>%
  inner_join(genes_map, by = "RSID") %>%
  filter(grepl(paste(whitelist_keywords, collapse="|"), EFO_Trait, ignore.case = TRUE)) %>%
  filter(!grepl("hernia|hearing|body mass|stature|intelligence|intelligence|hernia", EFO_Trait, ignore.case=TRUE)) %>%
  filter(P_value < 1e-6) %>%
  mutate(logP = -log10(P_value)) %>%
  mutate(logP = ifelse(logP > 100, 100, logP))

# Filter to top traits for a clean horizontal table
top_traits <- plot_data %>%
  group_by(EFO_Trait) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  head(25) %>%
  pull(EFO_Trait)

plot_data_sub <- plot_data %>%
  filter(EFO_Trait %in% top_traits)

# 3. Visualization
message("Generating Linear Pleiotropy Map...")

# Professional Bubbles
p <- ggplot(plot_data_sub, aes(x = Gene, y = reorder(EFO_Trait, logP))) +
  geom_point(aes(size = logP, fill = logP), shape = 21, color = "black", alpha = 0.8) +
  scale_fill_viridis_c(option = "magma", name = "-log10(P)") +
  scale_size_continuous(range = c(2, 10), name = "Association Significance") +
  theme_bw() +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.grid.major = element_line(color = "grey90", size = 0.1),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "right"
  ) +
  labs(
    title = "SLE Genetic Cross-Talk Portfolio",
    subtitle = "Shared susceptibility between SLE loci and other autoimmune/blood traits (GWAS Catalog)",
    x = "SLE High-Confidence Risk Genes",
    y = "Pleiotropic Trait (All labels horizontal and readable)"
  )

# Save
dir.create("figures", showWarnings = FALSE)
ggsave("figures/pleiotropy_map_linear.png", p, width = 12, height = 10, dpi = 300)
ggsave("figures/pleiotropy_map_linear.pdf", p, width = 12, height = 10)

message("Linear Pleiotropy Map saved to figures/pleiotropy_map_linear.png/pdf")
