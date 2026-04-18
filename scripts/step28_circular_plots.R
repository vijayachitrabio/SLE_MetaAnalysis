#!/usr/bin/env Rscript
# scripts/step28_circular_plots.R
# High-aesthetic Circular Manhattan and Pleiotropy Chord Diagram

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(circlize)
  library(RColorBrewer)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# ---------------------------------------------------------
# MODULE A: Circular Manhattan Plot
# ---------------------------------------------------------
message("=== Generating Circular Manhattan Plot ===")

# Load data (Downsample for performance, keeping discovery loci)
# We use top_merged_significant_loci from Step 21 as anchors
master <- fread("results/master_results_table.tsv")
discovery_loci <- master %>% filter(Final_Assessment == "HIGH CONFIDENCE")

# Load a representative subset of the full results
# If the 1.3GB file is too slow, we'll downsample significantly
results_file <- "results/discovery_meta_results.tsv"
if(file.exists(results_file)){
  # Read a subset for background points
  full_res <- fread(results_file, select = c("CHR", "BP", "P_meta", "RSID"))
  # Downsample background (keep 1% or ~20k SNPs)
  bg_data <- full_res[sample(.N, min(.N, 30000))]
  # Ensure discovery loci are included
  plot_data <- rbind(bg_data, full_res[RSID %in% discovery_loci$RSID]) %>% unique()
} else {
  # Fallback to master results if full results missing
  plot_data <- master %>% select(CHR, BP, P_meta, RSID)
}

plot_data <- plot_data %>%
  mutate(minuslog10p = -log10(P_meta)) %>%
  filter(!is.na(CHR))

# Prepare coordinates
nCHR <- length(unique(plot_data$CHR))
plot_data$CHR <- as.numeric(plot_data$CHR)
plot_data <- plot_data %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% 
  select(-chr_len) %>%
  left_join(plot_data, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

# Circular Manhattan using ggplot2 + coord_polar
p_circ <- ggplot(plot_data, aes(x=BPcum, y=minuslog10p)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.6, size=0.8) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype="dashed", alpha=0.5) +
  scale_color_manual(values = rep(c("#2c3e50", "#2980b9"), 22)) +
  coord_polar() +
  theme_void() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5, face="bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, size=10)
  ) +
  labs(title = "SLE Genetic Architecture", subtitle = "Circular Manhattan Mapping (European Ancestry)")

ggsave("figures/circular_manhattan.png", p_circ, width = 10, height = 10, dpi = 300)

# ---------------------------------------------------------
# MODULE B: Pleiotropy Chord Diagram
# ---------------------------------------------------------
message("=== Generating Pleiotropy Chord Diagram ===")

phewas <- fread("results/phewas_summary_refined.tsv")
top_loci <- master %>% filter(Final_Assessment == "HIGH CONFIDENCE") %>% select(RSID, Gene) %>% unique()

# Prepare link data
links <- phewas %>%
  inner_join(top_loci, by = "RSID") %>%
  filter(Category == "Immune-Mediated") %>%
  filter(P_value < 1e-6) %>%
  mutate(Gene = ifelse(Gene == "" | is.na(Gene), RSID, Gene)) %>%
  select(Gene, EFO_Trait) %>%
  group_by(Gene, EFO_Trait) %>%
  summarize(value = n()) %>%
  ungroup()

# Filter to top connections to keep it clean
top_links <- links %>%
  group_by(EFO_Trait) %>%
  mutate(total = sum(value)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  head(30)

# Setup colors
grid_col <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Set3"))
names(grid_col) <- unique(c(top_links$Gene, top_links$EFO_Trait))

png("figures/pleiotropy_chord_diagram.png", width = 1000, height = 1000, res = 150)
circos.par(gap.after = c(rep(2, length(unique(top_links$Gene))-1), 10, rep(2, length(unique(top_links$EFO_Trait))-1), 10))
chordDiagram(top_links, grid.col = grid_col, transparency = 0.4, 
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(top_links))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
}, bg.border = NA)
title("Genetic Pleiotropy: SLE Loci connectivity to Immune Traits")
dev.off()
circos.clear()

message("Done. Visuals saved to figures/.")
