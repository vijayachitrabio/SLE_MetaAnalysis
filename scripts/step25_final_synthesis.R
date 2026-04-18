#!/usr/bin/env Rscript
# scripts/step25_final_synthesis.R
# Final Synthesis of SLE Genetic Architecture
# Merging LAVA, COLOC, Meta-analysis, and Therapeutic Mapping

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("=========================================\n")
cat("Starting Final Results Synthesis\n")
cat("=========================================\n")

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# 1. Load Datasets
# Main Loci list (from audited LAVA run)
loci <- fread("results/final_lava_consolidated_loci.tsv")

# Colocalization Results
coloc <- fread("results/coloc_results_summary.tsv")

# Therapeutic Mapping
drug <- fread("results/therapeutic_mapping_summary.tsv")

# Novelty Details
novelty <- fread("results/precise_novelty_ensembl.tsv")

# 2. Integrate Colocalization (Effector Genes)
# Aggregate COLOC to one row per locus
if (nrow(coloc) > 0) {
  coloc_agg <- coloc %>%
    group_by(Locus) %>%
    summarise(
      Effector_Gene = paste(unique(Gene), collapse = ", "),
      Best_PP4 = max(PP4),
      Coloc_Tissue = paste(unique(Tissue), collapse = ", ")
    ) %>%
    ungroup()
  setDT(coloc_agg)
  loci <- merge(loci, coloc_agg, by.x = "RSID", by.y = "Locus", all.x = TRUE)
} else {
  loci[, `:=`(Effector_Gene = NA, Best_PP4 = NA, Coloc_Tissue = NA)]
}

# 3. Integrate Drug Targets
drug_subset <- drug[, .(RSID, Specific_Drug, Drug_Status, Drug_Class)]
loci <- merge(loci, drug_subset, by = "RSID", all.x = TRUE)

# 4. Integrate Precise Novelty labels
# In final_lava_consolidated_loci, we already have a 'Novelty' column, 
# but we can refine it if 'precise_novelty_ensembl' has better tags.
# Actually, let's just make sure we prioritize the 'HIGH CONFIDENCE' assessment.

# 5. Logical High-Confidence Summary
# Already defined in step23, let's just clean it up
loci[, Functional_Evidence := "LAVA (Genetic Correlation)"]
loci[!is.na(Best_PP4) & Best_PP4 > 0.8, Functional_Evidence := "LAVA + COLOC (eQTL)"]

# Sort by significance
loci <- loci[order(P_meta)]

# 6. Save Master Table
fwrite(loci, "results/master_results_table.tsv", sep="\t")

# 7. Print Final Publication Summary
cat("\n=== SLE Meta-Analysis: Final Architecture Summary ===\n")
cat(sprintf("Total Significant Loci: %d\n", nrow(loci)))
cat(sprintf("Novel Loci (Putative):   %d\n", sum(loci$Novelty == "Putative Novel")))
cat(sprintf("High-Confidence Loci:    %d\n", sum(loci$Final_Assessment == "HIGH CONFIDENCE")))
cat(sprintf("Loci with COLOC Evidence: %d\n", sum(!is.na(loci$Best_PP4))))

cat("\nTop High-Confidence Novel Loci:\n")
print(head(loci[Final_Assessment == "HIGH CONFIDENCE" & Novelty == "Putative Novel"], 5))

cat("\nSynthesis Complete! File saved to: results/master_results_table.tsv\n")
