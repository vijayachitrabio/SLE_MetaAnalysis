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
therap_rep <- therapeutic[Replicated == TRUE]

# Drug class colors - simplified labels
drug_colors <- c(
  "JAK Inhibitors" = "#E91E63",
  "Anti-IFNAR1 mAb" = "#9C27B0",
  "Complement Inhibitors" = "#3F51B5",
  "BAFF inhibitors" = "#009688",
  "TNF Pathway" = "#4CAF50",
  "Anti-CD20" = "#8BC34A",
  "Calcineurin Inhibitors" = "#CDDC39",
  "Anti-IFN" = "#FFEB3B",
  "TYK2 Inhibitors" = "#FFC107",
  "Anti-OX40L" = "#FF9800",
  "IVIG modulation" = "#795548",
  "Anti-IL-12/23" = "#607D8B",
  "Other" = "#9E9E9E"
)

# Simplify drug class names for plotting
therap_rep[, Drug_Class_simple := sub(" / .*", "", Drug_Class)]
therap_rep[, Drug_Class_simple := gsub("NF-kB.*", "TNF Pathway", Drug_Class_simple)]
therap_rep[, Drug_Class_simple := gsub("BAFF inhibitors.*", "BAFF inhibitors", Drug_Class_simple)]
therap_rep[, Drug_Class_simple := gsub("Anti-CD20 /.*", "Anti-CD20", Drug_Class_simple)]
therap_rep[, Drug_Class_simple := gsub("Anti-IFN /.*", "Anti-IFN", Drug_Class_simple)]
therap_rep[, Drug_Class_simple := gsub("IVIG /.*", "IVIG modulation", Drug_Class_simple)]

# Create a cleaner horizontal version
p_therap <- ggplot(therap_rep, aes(x = Drug_Class_simple, y = reorder(Lead_Gene, -Discovery_P))) +
  geom_point(aes(size = -log10(Discovery_P), color = Drug_Class_simple), alpha = 0.85, shape = 16) +
  scale_color_manual(values = drug_colors, name = "Drug Class") +
  scale_size_continuous(name = "-log10(P)", range = c(3, 10)) +
  labs(title = "Therapeutic Mapping: Drug Repurposing Potential",
       subtitle = "Replicated SLE loci with existing drug targets",
       x = "", y = "Target Gene") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(face = "italic", size = 10),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 1))

ggsave("figures/therapeutic_mapping.png", p_therap, width = 12, height = 10, dpi = 300, bg = "white")
ggsave("figures/therapeutic_mapping.pdf", p_therap, width = 12, height = 10)

# Create drug status summary
drug_summary <- therap_rep %>%
  group_by(Drug_Status) %>%
  summarise(n = n(), genes = paste(Lead_Gene, collapse = ", ")) %>%
  ungroup()

fwrite(drug_summary, "results/drug_status_summary.tsv", sep = "\t")

# ============================================================
# 3. TOP LOCI SUMMARY TABLE FIGURE
# ============================================================

loci <- fread("results/top_loci_summary_table.tsv")
loci_plot <- loci %>%
  filter(!is.na(Gene) & Gene != "TBD" & Gene != "Intergenic") %>%
  mutate(Gene_clean = sapply(strsplit(Gene, " / "), function(x) tail(x, 1))) %>%
  head(30)

p_table <- ggplot(loci_plot, aes(x = reorder(Gene_clean, P_disco), y = -log10(P_disco))) +
  geom_bar(aes(fill = Replicated), stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#4CAF50", "FALSE" = "#FF5722"), name = "Replicated") +
  coord_flip() +
  labs(title = "Top 30 SLE-Associated Loci",
       subtitle = "Sorted by discovery P-value",
       x = "", y = expression(-log[10](P[discovery]))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/top_loci_summary.png", p_table, width = 10, height = 12, dpi = 300, bg = "white")
ggsave("figures/top_loci_summary.pdf", p_table, width = 10, height = 12)

# ============================================================
# 4. NOVELTY CLASSIFICATION VISUALIZATION
# ============================================================

novelty_counts <- loci %>%
  filter(!is.na(Known_SLE)) %>%
  mutate(Novelty = ifelse(Known_SLE == TRUE, "Known SLE", "Novel")) %>%
  group_by(Novelty) %>%
  summarise(n = n())

p_novelty <- ggplot(novelty_counts, aes(x = Novelty, y = n, fill = Novelty)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_text(aes(label = n), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2", name = "Novelty Status") +
  labs(title = "Novelty Classification of SLE Loci",
       subtitle = "Known vs Novel genetic associations",
       x = "", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/novelty_classification.png", p_novelty, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("figures/novelty_classification.pdf", p_novelty, width = 6, height = 5)

cat("All downstream visualizations complete.\n")
cat("Output files in figures/:\n")
cat("- SA1_heterogeneity.png/pdf\n")
cat("- SA2_hla_distance.png\n")
cat("- SA3_replication.png/pdf\n")
cat("- therapeutic_mapping.png/pdf\n")
cat("- top_loci_summary.png/pdf\n")
cat("- novelty_classification.png/pdf\n")