#!/usr/bin/env Rscript
# plot_pathway_categories.R (v4)
# aesthetically refined categorized pathway plot for publication.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

# --- Configuration ---
DIR_IN  <- "results/pathway"
DIR_OUT <- "figures"
FILE_FGSEA <- file.path(DIR_IN, "fgsea_significant_all.tsv")
OUTPUT_PNG <- file.path(DIR_OUT, "Fig_2_pathway_functional_v4.png")
OUTPUT_PDF <- file.path(DIR_OUT, "Fig_2_pathway_functional_v4.pdf")

# --- Load and Categorize Data ---
cat("Loading and processing pathway data...\n")
fgsea_dt <- fread(FILE_FGSEA)

categorize_term <- function(term) {
  term <- toupper(term)
  if (grepl("COMPLEMENT|OPSONIZATION|PHAGOCYTOSIS", term)) return("Complement & Phagocytosis")
  if (grepl("MHC|ANTIGEN|ADAPTIVE|HLA|TAP_BINDING", term)) return("Antigen Presentation & HLA")
  if (grepl("NATURAL_KILLER|NK_CELL", term)) return("NK Cell Biology")
  if (grepl("INTERFERON", term)) return("Interferon Signaling")
  if (grepl("APOPTOTIC|NECROPTOTIC|CELL_CLEARANCE", term)) return("Cell Death & Clearance")
  if (grepl("URATE|SLC22|SLC17", term)) return("Transport & Metabolism")
  return("Other Immune / Structural")
}

# Process results
plot_dt <- fgsea_dt %>%
  mutate(
    category = sapply(pathway, categorize_term),
    # Cleaner pathway labels
    short_name = gsub("^(GOBP|REACTOME|HALLMARK)_", "", pathway),
    short_name = gsub("_", " ", short_name),
    short_name = str_to_sentence(short_name),
    short_name = case_when(
      short_name == "Adaptive immune response" ~ "Adaptive Immune Response",
      short_name == "Interferon gamma response" ~ "IFN-gamma Response",
      TRUE ~ short_name
    ),
    log10p = -log10(padj)
  ) %>%
  # Filter to focus on key biological story
  filter(category != "Other Immune / Structural" | short_name == "Adaptive Immune Response")

# Set Category and Term factor order
cat_order <- c("Complement & Phagocytosis", "Antigen Presentation & HLA", 
               "NK Cell Biology", "Interferon Signaling", 
               "Cell Death & Clearance", "Transport & Metabolism")
plot_dt$category <- factor(plot_dt$category, levels = cat_order)

plot_dt <- plot_dt %>%
  arrange(category, NES) %>%
  mutate(short_name = factor(short_name, levels = unique(short_name)))

# --- Generate Visualization ---
cat("Generating refined categorized bubble plot...\n")

p <- ggplot(plot_dt, aes(x = NES, y = short_name)) +
  # Background Shading for Categories (using geom_rect)
  # To do this cleanly with facets, we'll use facet headers or a custom strip
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  # Points: Size = Genes, Fill = Significance
  geom_point(aes(size = size, fill = log10p), shape = 21, color = "black", stroke = 0.2) +
  # Color Scale: Dark Red for high significance
  scale_fill_gradient(low = "#fee5d9", high = "#99000d", 
                      name = expression(-log[10](P[adj]))) +
  # Size Scale: Number of genes
  scale_size_continuous(range = c(4, 12), name = "Number of Genes") +
  # Faceting with category labels on the right
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Functionally Categorized Pathway Enrichment",
       subtitle = "SLE European Meta-Analysis | bubble size = gene count | bubble color = significance",
       x = "Normalised Enrichment Score (NES)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    # Professional Theme Elements
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey98"),
    
    # Category Headers (Strip)
    strip.background = element_rect(fill = "#fdfdfd", color = "grey90"),
    strip.text.y = element_text(angle = 0, face = "bold", size = 10, hjust = 0),
    
    # Axis & Labels
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    
    # Legend at Bottom
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    
    # Overall plot
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b = 15)),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = 10, barheight = 1),
    size = guide_legend(order = 2)
  )

# Save
cat("Saving high-resolution plot to PNG and PDF formats...\n")
ggsave(OUTPUT_PNG, p, width = 12, height = 10, dpi = 350, bg = "white")
ggsave(OUTPUT_PDF, p, width = 12, height = 10, device = "pdf", bg = "white")

cat("Done.\n")
