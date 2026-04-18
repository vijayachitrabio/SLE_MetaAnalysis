#!/usr/bin/env Rscript
# scripts/step11_locus_plots.R
# Generate regional association plots for Top 5 SLE Discovery Loci.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

DIR_RESULTS <- "results"
FILE_META    <- file.path(DIR_RESULTS, "discovery_meta_results.tsv")
DIR_FIGS     <- "figures/locus_plots"
dir.create(DIR_FIGS, showWarnings = FALSE, recursive = TRUE)

# 1. Load Audited High-Confidence Loci to Plot (from Step 21 summary)
top_loci <- fread("results/top_loci_summary_table.tsv")
loci_to_plot <- top_loci %>%
  filter(Replicated == TRUE) %>%
  select(Gene, RSID, CHR, BP) %>%
  arrange(desc(Gene != "")) %>% # Prioritize loci with gene names
  head(15)

# Clean characters that might break filenames
loci_to_plot$Gene <- gsub("/", "-", loci_to_plot$Gene)
loci_to_plot$Gene[loci_to_plot$Gene == ""] <- paste0("Locus_", loci_to_plot$RSID[loci_to_plot$Gene == ""])

cat("=== Generating Regional Locus Plots ===\n")

# 2. Function to plot a region
plot_locus <- function(gene, rsid, chr, bp) {
  window <- 500000 # 500kb each side (1MB total)
  start_loc <- bp - window
  end_loc   <- bp + window
  
  cat(sprintf("  Plotting %s (%s) at Chr%d:%d-%d...\n", gene, rsid, chr, start_loc, end_loc))
  
  # Read meta-results for this region only
  # Efficiently read using awk or fread filtering if possible, 
  # but since it's 1.3GB we'll just read and filter for now as it's manageable for a 1MB window if indexed.
  # Here we'll read the whole file but only once if we were smart, 
  # but for simplicity we'll read it once and iterate.
  
  region_data <- fread(FILE_META) %>%
    filter(CHR == chr, BP >= start_loc, BP <= end_loc)
  
  if(nrow(region_data) == 0) {
    cat(sprintf("    Warning: No data found for %s region.\n", gene))
    return(NULL)
  }
  
  region_data <- region_data %>%
    mutate(is_lead = (RSID == rsid)) %>%
    mutate(minuslog10p = -log10(P_meta))
  
  p <- ggplot(region_data, aes(x = BP / 1e6, y = minuslog10p)) +
    geom_point(aes(color = is_lead, size = is_lead), alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#e74c3c")) +
    scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 4)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = subset(region_data, is_lead), aes(label = RSID), 
                    nudge_y = 2, fontface = "bold", color = "#e74c3c") +
    theme_minimal() +
    labs(title = paste("Regional Association:", gene),
         subtitle = paste("Lead SNP:", rsid, "| Chromosome", chr),
         x = paste("Position on Chromosome", chr, "(Mb)"),
         y = expression(-log[10](P[meta]))) +
    theme(legend.position = "none")
  
  ggsave(file.path(DIR_FIGS, paste0("locus_", gene, ".png")), p, width = 8, height = 5, dpi = 300)
}

# 3. Process each locus
# Reading 1.3GB once and then subsetting is much faster
cat("Reading meta-analysis results into memory for plotting...\n")
all_meta <- fread(FILE_META)

for (i in 1:nrow(loci_to_plot)) {
  gene <- loci_to_plot$Gene[i]
  rsid <- loci_to_plot$RSID[i]
  chr  <- loci_to_plot$CHR[i]
  bp   <- loci_to_plot$BP[i]
  
  # Subsetting logic
  window <- 500000
  start_loc <- bp - window
  end_loc   <- bp + window
  
  region_data <- all_meta[CHR == chr & BP >= start_loc & BP <= end_loc]
  
  if(nrow(region_data) == 0) next
  
  region_data <- region_data %>%
    mutate(is_lead = (RSID == rsid)) %>%
    mutate(minuslog10p = -log10(P_meta))
  
  p <- ggplot(region_data, aes(x = BP / 1e6, y = minuslog10p)) +
    geom_point(aes(color = is_lead, size = is_lead, alpha = is_lead)) +
    scale_color_manual(values = c("FALSE" = "#bdc3c7", "TRUE" = "#e74c3c")) +
    scale_size_manual(values = c("FALSE" = 1.0, "TRUE" = 3.5)) +
    scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1.0)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#34495e") +
    geom_text_repel(data = subset(region_data, is_lead), aes(label = RSID), 
                    nudge_y = 3, fontface = "bold", color = "#e74c3c", segment.color = "#e74c3c") +
    theme_classic() +
    labs(title = paste("Association Locus:", gene),
         subtitle = paste("Lead SNP:", rsid, "| Build GRCh38"),
         x = paste("Chromosome", chr, "Position (Mb)"),
         y = expression(-log[10](P))) +
    theme(legend.position = "none", 
          plot.title = element_text(face="bold", size=16),
          axis.title = element_text(face="bold"))
  
  ggsave(file.path(DIR_FIGS, paste0("locus_", gene, ".png")), p, width = 8, height = 5, dpi = 300)
}

cat("Locus plots generated successfully in figures/locus_plots/.\n")
