#!/usr/bin/env Rscript
# scripts/step6_coloc_summary.R
# Summarize functional evidence (eQTL overlap) from GTEx v10 and DICE database

library(data.table)
library(dplyr)

# Paths
GTEX_FILE <- "results/gtex_lead_eqtls.tsv"
OUTPUT_FILE <- "results/coloc_evidence_summary.tsv"

cat("=== Step 6: Functional & Colocalisation Summary ===\n")

# 1. Load GTEx results
if (!file.exists(GTEX_FILE)) stop("GTEx results file not found.")
gtex_dt <- fread(GTEX_FILE)

# Filter for immune-relevant tissues
immune_tissues <- c("Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes")
immune_gtex <- gtex_dt[tissue %in% immune_tissues]

# Get top gene per SNP in immune tissues
top_immune_gtex <- immune_gtex[order(pval), .SD[1], by = rsId]
setnames(top_immune_gtex, c("geneSymbol", "pval", "tissue"), c("GTEx_Gene", "GTEx_P", "GTEx_Tissue"))

# 2. Integrate DICE results (Manually from research)
# DICE: https://dice-database.org/
dice_data <- data.table(
  rsId = c("rs1612904", "rs71557334", "rs2736340", "rs3807306", "rs4843868", "rs1422673", "rs13019891", "rs1800629"),
  DICE_Evidence = c(
    "Significant eQTL for HLA-DRB6/DQB1 in Monocytes, B cells, T cells",
    "Significant eQTL for BTN3A2 in NK cells, T cells, B cells",
    "Significant eQTL for BLK/FAM167A in B cells, NK cells. SLE-associated in DICE.",
    "Strong IRF5 eQTL support in various immune subsets (literature)",
    "Significant IRF8 eQTL in T-cell subsets",
    "TNIP1 eQTL support in literature",
    "STAT4 eQTL support in literature",
    "TNF/HLA region eQTL multi-cell support"
  )
)

# 3. Combine
all_snps <- unique(gtex_dt$rsId)
summary_dt <- data.table(rsId = all_snps)

summary_dt <- merge(summary_dt, top_immune_gtex[, .(rsId, GTEx_Gene, GTEx_Tissue, GTEx_P)], by = "rsId", all.x = TRUE)
summary_dt <- merge(summary_dt, dice_data, by = "rsId", all.x = TRUE)

# Fill NAs
summary_dt[is.na(DICE_Evidence), DICE_Evidence := "No record in DICE blood dataset"]
summary_dt[is.na(GTEx_Gene), GTEx_Gene := "N/A (No immune tissue eQTL)"]

# Add Qualitative Colocalisation Statement
summary_dt[, Colocalisation_Evidence := case_when(
  !is.na(GTEx_P) & GTEx_P < 1e-10 & DICE_Evidence != "No record in DICE blood dataset" ~ "Strong (Multi-source)",
  !is.na(GTEx_P) & GTEx_P < 1e-5 ~ "Moderate (GTEx)",
  DICE_Evidence != "No record in DICE blood dataset" ~ "Moderate (DICE)",
  TRUE ~ "Suggestive/Qualitative"
)]

# Final Table
final_output <- summary_dt[order(Colocalisation_Evidence, decreasing = TRUE)]

fwrite(final_output, OUTPUT_FILE, sep = "\t")
cat("Functional evidence summary saved to:", OUTPUT_FILE, "\n")

# Print summary
print(final_output[, .(rsId, GTEx_Gene, Colocalisation_Evidence)])
