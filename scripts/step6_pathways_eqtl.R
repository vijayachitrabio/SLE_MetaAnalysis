library(data.table)
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(ggpubr)
library(stringr)

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

message("Loading full audited loci mapping...")
master <- fread("results/master_results_table.tsv")

# Extract unique genes (Effector_Gene prioritized, then Gene)
genes <- master %>%
  mutate(Gene_Final = coalesce(na_if(Effector_Gene, ""), na_if(Gene, ""))) %>%
  filter(!is.na(Gene_Final) & Gene_Final != "" & Gene_Final != "TBD" & Gene_Final != "Intergenic") %>%
  pull(Gene_Final) %>%
  unique()

message(paste("Found", length(genes), "unique mapped genes for enrichment."))

message("Running gProfiler enrichment...")
gostres <- gost(query = genes, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO:BP", "REAC", "KEGG"))

if (is.null(gostres)) {
  message("No significant pathways found.")
} else {
  message("Processing and filtering enrichment results...")
  
  enrich_df <- gostres$result %>% 
    select(source, term_id, term_name, p_value, intersection_size, query_size) %>%
    mutate(
      log_pval = -log10(p_value),
      overlap_ratio = intersection_size / query_size * 100,
      # Clean term names
      term_clean = str_remove_all(term_name, "^GO:.*? "),
      term_clean = ifelse(nchar(term_clean) > 60, str_trunc(term_clean, 57), term_clean)
    ) %>%
    arrange(p_value)
  
  # Filter for non-redundant, most relevant terms
  # Prioritize GO:BP and REAC with good intersection sizes
  enrich_filtered <- enrich_df %>%
    filter(source %in% c("GO:BP", "REAC")) %>%
    filter(intersection_size >= 3) %>%
    group_by(source) %>%
    slice_head(n = 15) %>%
    ungroup() %>%
    arrange(p_value) %>%
    head(25)
  
  # Save full results
  fwrite(enrich_df, "results/pathway_enrichment_results.tsv", sep="\t")
  
  # Create enhanced visualization
  message("Creating pathway enrichment visualization...")
  
  # Color palette for sources
  source_colors <- c("GO:BP" = "#4CAF50", "REAC" = "#2196F3", "KEGG" = "#FF9800")
  
  # Dot plot - top pathways
  p1 <- ggplot(enrich_filtered, aes(x = reorder(term_clean, -p_value), y = log_pval)) +
    geom_point(aes(size = intersection_size, color = source), alpha = 0.7) +
    scale_color_manual(values = source_colors, name = "Source") +
    scale_size_continuous(name = "Genes", range = c(3, 10)) +
    coord_flip() +
    labs(title = "Pathway Enrichment Analysis: SLE-associated Loci",
         subtitle = "Top 25 significant pathways (FDR < 0.05)",
         x = "", y = expression(-log[10](P))) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")
  
  # Bar plot with -log10 p-value
  p2 <- ggplot(enrich_filtered %>% head(20), 
               aes(x = reorder(term_clean, log_pval), y = log_pval, fill = source)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = source_colors) +
    coord_flip() +
    labs(title = "Top 20 Enriched Pathways",
         subtitle = "SLE-associated genes",
         x = "", y = expression(-log[10](P))) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Combine plots
  p_combined <- ggarrange(p1, p2, ncol = 2, labels = c("A", "B"), 
                          font.label = list(size = 12, face = "bold"))
  
  # Save outputs
  ggsave("figures/pathway_enrichment.png", p1, width = 12, height = 10, dpi = 300, bg = "white")
  ggsave("figures/pathway_enrichment_barplot.png", p2, width = 10, height = 8, dpi = 300, bg = "white")
  ggsave("figures/pathway_enrichment_combined.pdf", p_combined, width = 16, height = 10)
  ggsave("figures/pathway_enrichment.pdf", p1, width = 12, height = 10)
  
  # Save filtered table
  fwrite(enrich_filtered, "results/pathway_enrichment_filtered.tsv", sep = "\t")
  
  message(paste("Saved", nrow(enrich_filtered), "filtered pathways to results/pathway_enrichment_filtered.tsv"))
  message("Done with step 6.")
}
