# scripts/step23_consolidate_lava.R
library(data.table)
library(dplyr)

# 1. Load Data
mapping <- fread("results/top_loci_lava_mapping.tsv")
lava <- fread("results/lava_sle_results.csv")

# 2. Merge
# We merge by LOC
merged <- merge(mapping, lava, by.x = "LAVA_LOC", by.y = "LOC", all.x = TRUE)

# 3. Categorize
# Threshold for LAVA confirmation: P < 0.05 (nominal) or P < 0.05 / 50 (Bonferroni for interesting blocks)
# Actually, LAVA P is the significance of the genetic correlation.
merged[, LAVA_Status := "Filtered (Spurious/No Local Signal)"]
merged[!is.na(p) & p < 0.05, LAVA_Status := "Confirmed (Nominal)"]
merged[!is.na(p) & p < (0.05/50), LAVA_Status := "Confirmed (Bonferroni)"]

# Handle loci where LAVA couldn't run (missing p)
merged[is.na(p), LAVA_Status := "Insufficient Local Signal"]

# 4. Find Redundant Loci
loc_counts <- merged[, .N, by = LAVA_LOC]
redundant_locs <- loc_counts[N > 1 & !is.na(LAVA_LOC), LAVA_LOC]
merged[, Redundancy := "Independent"]
merged[LAVA_LOC %in% redundant_locs, Redundancy := "Shared LAVA Block"]

# 5. Final Combined Status
merged[, Final_Assessment := ifelse(LAVA_Status %like% "Confirmed", "HIGH CONFIDENCE", "LOW CONFIDENCE/SKEPTICAL")]

# 6. Save Updated Table
final_output <- merged[, .(
  RSID, CHR, BP, P_meta, BETA, 
  LAVA_LOC, LAVA_rg = rho, LAVA_P = p,
  LAVA_Status, Redundancy, Final_Assessment,
  Gene, Region, Novelty
)]

fwrite(final_output, "results/final_lava_consolidated_loci.tsv", sep="\t")

# Summary Stats
cat("\n=== LAVA Consolidation Summary ===\n")
cat("Total Candidate Loci:", nrow(final_output), "\n")
cat("High Confidence (Confirmed by LAVA):", sum(final_output$Final_Assessment == "HIGH CONFIDENCE"), "\n")
cat("Low Confidence (Filtered/Weak):", sum(final_output$Final_Assessment == "LOW CONFIDENCE/SKEPTICAL"), "\n")
cat("\nTop 10 High Confidence Loci:\n")
print(head(final_output[Final_Assessment == "HIGH CONFIDENCE"][order(P_meta)], 10))
