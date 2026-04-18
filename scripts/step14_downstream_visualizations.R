#!/usr/bin/env Rscript
# scripts/step14_downstream_visualizations.R
# Additional downstream visualizations: Sensitivity analyses, Therapeutic mapping, Summary tables

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(stringr)
  library(igraph)
  library(ggrepel)
  library(RColorBrewer)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

cat("=== Creating Downstream Visualizations ===\n")

# ============================================================
# 1. SENSITIVITY ANALYSIS PLOTS
# ============================================================

# SA1: Random vs Fixed Effects (Heterogeneity)
sa1 <- fread("results/sensitivity/SA1_random_vs_fixed.tsv")
sa1[, Note_clean := factor(Note, levels = c("Consistent", "No between-study variance"))]

p_sa1 <- ggplot(sa1, aes(x = reorder(Gene, -P_fixed), y = I2)) +
  geom_bar(fill = "steelblue", alpha = 0.8, stat = "identity") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 75, linetype = "dotted", color = "darkred", linewidth = 1) +
  coord_flip() +
  labs(title = "SA1: Heterogeneity (I²) by Locus",
       subtitle = "Fixed-effects meta-analysis",
       x = "", y = "I² (%)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  annotate("text", x = 5, y = 52, label = "Moderate (50%)", color = "red", size = 3) +
  annotate("text", x = 5, y = 77, label = "High (75%)", color = "darkred", size = 3)

ggsave("figures/SA1_heterogeneity.png", p_sa1, width = 10, height = 12, dpi = 300, bg = "white")
ggsave("figures/SA1_heterogeneity.pdf", p_sa1, width = 10, height = 12)

# SA3: Replication Results
sa3 <- fread("results/sensitivity/SA3_replication_power.tsv")
sa3[, Conclusion_clean := factor(Conclusion, levels = c("Replicated", "Adequately powered - non-replicating", "Underpowered"))]

p_sa3 <- ggplot(sa3, aes(x = reorder(Gene, Power_pct), y = Power_pct, fill = Conclusion_clean)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("Replicated" = "#2E7D32", 
                               "Adequately powered - non-replicating" = "#FF9800",
                               "Underpowered" = "#9E9E9E"),
                   name = "Replication Status") +
  coord_flip() +
  labs(title = "SA3: Replication Analysis",
       subtitle = "Power and replication status in Spanish cohort",
       x = "", y = "Statistical Power (%)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/SA3_replication.png", p_sa3, width = 10, height = 14, dpi = 300, bg = "white")
ggsave("figures/SA3_replication.pdf", p_sa3, width = 10, height = 14)

# SA2: HLA Pairwise Distance (Network)
sa2 <- fread("results/sensitivity/SA2_HLA_pairwise.tsv")
hla_genes <- unique(c(sa2$Gene_A, sa2$Gene_B))
hla_genes <- hla_genes[grepl("HLA", hla_genes)]

# Create simple distance-based visualization
hla_dist <- sa2 %>% 
  filter(grepl("HLA", Gene_A) & grepl("HLA", Gene_B)) %>%
  mutate(Gene_A = sub(".* / ", "", Gene_A),
         Gene_B = sub(".* / ", "", Gene_B)) %>%
  select(Gene_A, Gene_B, Distance_kb)

p_sa2 <- ggplot(hla_dist, aes(x = reorder(Gene_A, -Distance_kb), y = Gene_B)) +
  geom_tile(aes(fill = Distance_kb), color = "white") +
  scale_fill_viridis_c(name = "Distance (kb)", direction = -1) +
  labs(title = "SA2: HLA Region Pairwise LD",
       subtitle = "Pairwise distance between lead SNPs",
       x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/SA2_hla_distance.png", p_sa2, width = 8, height = 6, dpi = 300, bg = "white")

cat("Sensitivity analysis plots saved\n")

# ============================================================
# 2. THERAPEUTIC MAPPING VISUALIZATION
# ============================================================

therapeutic <- fread("results/therapeutic_mapping_summary.tsv", nrows = 26)
therapeutic <- therapeutic[RSID != "" & !is.na(RSID)]

# Filter to replicated genes for cleaner visualization
# Filter to only High-Confidence Loci from the final summary
hc_loci <- fread("results/top_loci_summary_table.tsv")$RSID
therap_hc <- therapeutic[RSID %in% hc_loci]

# Custom professional color palette for Drug Status
status_colors <- c(
  "FDA-approved" = "#2E7D32",         # Royal Green
  "FDA-approved (RA/SLE)" = "#2E7D32", 
  "Approved" = "#2E7D32",
  "Approved / Phase III" = "#1565C0", # Deep Blue
  "Phase III" = "#1565C0",
  "Phase II" = "#F9A825",            # Amber/Orange
  "Phase II SLE" = "#F9A825",
  "Preclinical" = "#757575"           # Grey
)

# Clean up drug status for plotting
therap_hc[, Status_Simple := factor(Drug_Status, levels = c("FDA-approved", "FDA-approved (RA/SLE)", "Approved", "Approved / Phase III", "Phase III", "Phase II", "Phase II SLE", "Preclinical"))]

# Create the premium visualization
p_therap <- ggplot(therap_hc, aes(x = Drug_Class, y = reorder(Lead_Gene, -Discovery_P))) +
  # Background panel styling
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F5F5F5", alpha = 0.2) +
  # Bubble Plot
  geom_point(aes(size = -log10(Discovery_P), color = Drug_Status, shape = Drug_Status), alpha = 0.8) +
  # Labels with ggrepel
  geom_text_repel(aes(label = Specific_Drug), size = 3, fontface = "bold", box.padding = 0.5, max.overlaps = 15) +
  # Scales
  scale_color_manual(values = status_colors, name = "Clinical Status") +
  scale_shape_manual(values = c(16, 16, 16, 17, 17, 18, 18, 15), name = "Clinical Status") +
  scale_size_continuous(name = expression(-log[10](P[discovery])), range = c(4, 12)) +
  # Labels and Theme
  labs(title = "Figure 5: Therapeutic Landscape of SLE High-Confidence Loci",
       subtitle = "Mapping of validated genetic associations to existing and emerging drug targets",
       x = "Drug Category", y = "Target Effector Gene") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(face = "italic", size = 11),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, color = "#1A237E"),
    plot.subtitle = element_text(size = 11, color = "#455A64"),
    panel.grid.major = element_line(color = "white"),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "#FAFAFA", color = "grey90")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
  )

ggsave("figures/therapeutic_mapping.png", p_therap, width = 12, height = 10, dpi = 300, bg = "white")
ggsave("figures/therapeutic_mapping.pdf", p_therap, width = 12, height = 10)

# Create drug status summary
drug_summary <- therap_hc %>%
  group_by(Drug_Status) %>%
  summarise(n = n(), genes = paste(Lead_Gene, collapse = ", ")) %>%
  ungroup()

fwrite(drug_summary, "results/drug_status_summary.tsv", sep = "\t")

# ============================================================
# 3. TOP LOCI SUMMARY TABLE FIGURE
# ============================================================

loci <- fread("results/top_loci_summary_table.tsv")
loci_plot <- loci %>%
  filter(!is.na(Gene) & Gene != "TBD" & Gene != "Intergenic" & Gene != "") %>%
  mutate(Gene_clean = sapply(strsplit(Gene, " / "), function(x) tail(x, 1))) %>%
  head(30)

p_table <- ggplot(loci_plot, aes(x = reorder(Gene_clean, P_meta), y = -log10(P_meta))) +
  geom_bar(aes(fill = as.character(Replicated)), stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#4CAF50", "FALSE" = "#FF5722"), name = "Replicated") +
  coord_flip() +
  labs(title = "Top 30 SLE-Associated Loci",
       subtitle = "Sorted by discovery meta P-value",
       x = "", y = expression(-log[10](P[meta]))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/top_loci_summary.png", p_table, width = 10, height = 12, dpi = 300, bg = "white")
ggsave("figures/top_loci_summary.pdf", p_table, width = 10, height = 12)

# = [Removed Section 4: Novelty Classification] = #
cat("All downstream visualizations complete.\n")
cat("Output files in figures/:\n")
cat("- SA1_heterogeneity.png/pdf\n")
cat("- SA2_hla_distance.png\n")
cat("- SA3_replication.png/pdf\n")
cat("- therapeutic_mapping.png/pdf\n")
cat("- top_loci_summary.png/pdf\n")