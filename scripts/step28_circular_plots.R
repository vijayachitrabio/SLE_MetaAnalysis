#!/usr/bin/env Rscript
# scripts/step28_circular_plots.R
# Refined Pleiotropy Chord Diagram (No title, high-res PDF, optimized labels)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(circlize)
  library(RColorBrewer)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# 1. Load Data
phewas <- fread("results/phewas_summary_refined.tsv")
master <- fread("results/master_results_table.tsv")
top_loci <- master %>% 
  filter(Final_Assessment == "HIGH CONFIDENCE") %>% 
  select(RSID, Gene) %>% 
  unique()

# Prepare link data
links <- phewas %>%
  inner_join(top_loci, by = "RSID") %>%
  filter(Category == "Immune-Mediated") %>%
  filter(P_value < 5e-5) %>% # Slightly relaxed for better visual density
  mutate(Gene = ifelse(Gene == "" | is.na(Gene), RSID, Gene)) %>%
  select(Gene, EFO_Trait) %>%
  group_by(Gene, EFO_Trait) %>%
  summarize(value = n()) %>%
  ungroup()

# Filter to top connections to ensure readability
# Keep top 12 genes and top 12 traits for a clean circle
top_genes <- links %>% group_by(Gene) %>% summarize(v = sum(value)) %>% arrange(desc(v)) %>% head(12) %>% pull(Gene)
top_traits <- links %>% group_by(EFO_Trait) %>% summarize(v = sum(value)) %>% arrange(desc(v)) %>% head(12) %>% pull(EFO_Trait)

clean_links <- links %>%
  filter(Gene %in% top_genes & EFO_Trait %in% top_traits)

# 2. Setup Plotting Parameters
# Custom color palette
grid_col <- c(brewer.pal(length(unique(clean_links$Gene)), "Set1"), 
              brewer.pal(length(unique(clean_links$EFO_Trait)), "Set3"))
names(grid_col) <- unique(c(clean_links$Gene, clean_links$EFO_Trait))

# 3. Generating High-Resolution Plots (PNG and PDF)
message("Generating Refined Pleiotropy Chord Diagram...")

# Function to draw the diagram (Repeated for PNG and PDF)
draw_chord <- function() {
  circos.clear()
  # Adjust margins to ensure labels aren't cut off
  circos.par(
    gap.after = c(rep(2, length(unique(clean_links$Gene))-1), 10, rep(2, length(unique(clean_links$EFO_Trait))-1), 10),
    canvas.xlim = c(-1.3, 1.3),
    canvas.ylim = c(-1.3, 1.3)
  )
  
  chordDiagram(
    clean_links, 
    grid.col = grid_col, 
    transparency = 0.35, 
    annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.05)
  )
  
  # Custom labels with optimized text size and distancing
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ylim[1] + 1.2, 
      CELL_META$sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE, 
      adj = c(0, 0.5), 
      cex = 0.65, # Adjusted for uniform readability
      font = 2     # Bold labels
    )
  }, bg.border = NA)
}

# Save PNG (300 DPI)
png("figures/pleiotropy_chord_diagram.png", width = 1800, height = 1800, res = 300)
draw_chord()
dev.off()

# Save PDF (Vector format)
pdf("figures/pleiotropy_chord_diagram.pdf", width = 10, height = 10)
draw_chord()
dev.off()

circos.clear()
message("Done. Refined Pleiotropy Chord Diagram saved to figures/ (PNG and PDF).")
