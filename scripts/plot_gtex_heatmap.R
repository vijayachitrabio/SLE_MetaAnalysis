#!/usr/bin/env Rscript
# plot_gtex_heatmap.R (v2)
# Generate a high-quality multi-tissue eQTL heatmap for lead SLE variants.

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# --- Configuration ---
INPUT_FILE <- "results/gtex_lead_eqtls.tsv"
OUTPUT_PNG <- "figures/Fig_4_GTEx_heatmap.png"
OUTPUT_PDF <- "figures/Fig_4_GTEx_heatmap.pdf"

# --- Load and Filter Data ---
cat("Loading GTEx data...\n")
dt <- fread(INPUT_FILE)

# Define tissues of interest for the plot (focused view)
tissues_of_interest <- c(
  "Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes", 
  "Lung", "Skin_Not_Sun_Exposed_Suprapubic", "Nerve_Tibial",
  "Adipose_Subcutaneous", "Artery_Tibial", "Muscle_Skeletal", "Thyroid"
)

# Define primary SLE target genes to highlight (one per SNP for clarity)
plot_ready <- dt %>%
  filter(tissue %in% tissues_of_interest) %>%
  group_by(rsId) %>%
  # Priority: Canonical gene list then lowest p-value
  mutate(is_canonical = geneSymbol %in% c("STAT4", "IRF5", "TNF", "BLK", "XKR6", "C4A", "HLA-DRB1", "PAX8", "IL1A", "IRF8", "TNIP1", "ITGAX")) %>%
  arrange(rsId, desc(is_canonical), pval) %>%
  filter(geneSymbol == first(geneSymbol)) %>% 
  ungroup()

# Create clean labels
plot_ready <- plot_ready %>%
  mutate(row_label = paste0(rsId, " (", geneSymbol, ")")) %>%
  mutate(clean_tissue = str_replace_all(tissue, "_", " ") %>% str_wrap(25)) %>%
  # Significance logic
  mutate(stars = case_when(
    pval < 1e-10 ~ "***",
    pval < 1e-5  ~ "**",
    pval < 1e-3  ~ "*",
    TRUE ~ ""
  ))

# --- Plotting ---
cat("Generating premium heatmap...\n")

# Color palette: Blue (negative/downreg) -> White (neutral) -> Orange (positive/upreg)
# Using ColorBrewer RdBu inspired hex codes
color_low    <- "#0571b0" # Deep Blue
color_mid    <- "#f7f7f7" # Off-White
color_high   <- "#e66101" # Deep Orange

p <- ggplot(plot_ready, aes(x = clean_tissue, y = reorder(row_label, nes), fill = nes)) +
  # Use linewidth for separation
  geom_tile(color = "white", linewidth = 0.8) +
  # Significance stars: Bold and larger
  geom_text(aes(label = stars), color = "black", fontface = "bold", size = 6, vjust = 0.75) +
  # Diverging color scale
  scale_fill_gradient2(
    low = color_low, 
    mid = color_mid, 
    high = color_high, 
    midpoint = 0, 
    name = "Directional Effect\n(Normalized Effect Size)",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15, barheight = 1)
  ) +
  theme_minimal(base_family = "sans") +
  labs(
    title = "Tissue-Specific Expression Effects of SLE Risk Alleles",
    subtitle = "Directional impact on target gene expression across GTEx v10 tissues",
    x = NULL, 
    y = "Lead GWAS Variant (Target Gene)",
    caption = "Significance: *** P < 1e-10, ** P < 1e-5, * P < 1e-3\nBlue: Risk allele decreases expression; Orange: Risk allele increases expression."
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 20, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, color = "grey30", margin = margin(b = 20)),
    plot.caption = element_text(size = 10, color = "grey40", hjust = 0, margin = margin(t = 15)),
    plot.margin = margin(20, 20, 20, 20)
  )

# --- Save Outputs ---
cat("Saving heatmap to PNG and PDF formats...\n")
ggsave(OUTPUT_PNG, p, width = 13, height = 10, dpi = 350, bg = "white")
ggsave(OUTPUT_PDF, p, width = 13, height = 10, device = "pdf", bg = "white")

cat("Success.\n")
