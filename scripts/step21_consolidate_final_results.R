#!/usr/bin/env Rscript
# scripts/step21_consolidate_final_results.R
# Final consolidation of novelty, replication, and annotation

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

cat("=== Finalizing SLE Meta-Analysis Results ===\n")

# 1. Load the precise clumping results (from step20)
results <- fread("results/precise_novelty_ensembl.tsv")

# 2. Comprehensive Known Loci List (from step16)
known_loci <- data.table(
  RSID = c("rs2476601", "rs1145879", "rs12569394", "rs10954293", "rs2000999",
           "rs231735", "rs4963128", "rs11934158", "rs1320333", "rs2070190",
           "rs7573219", "rs11792679", "rs1079026", "rs3821236", "rs2010963",
           "rs10506418", "rs10181617", "rs10903045", "rs12737182", "rs1150754",
           "rs7214628", "rs2072114", "rs9380229", "rs9888739", "rs3027898",
           "rs2040406", "rs11078903", "rs11078902", "rs3795877", "rs2736340",
           "rs3099844", "rs9271546", "rs9271538", "rs3768797", "rs1150753",
           "rs10516526", "rs2662421", "rs11235604", "rs10806425", "rs2295420",
           "rs3127214", "rs2233945", "rs3134869", "rs3097673", "rs3097674",
           "rs9366772", "rs6931746", "rs2241769", "rs10800593", "rs10138183",
           "rs4853458", "rs35000415", "rs34572943", "rs10036748", "rs58721818"),
  Known_Status = "Known"
)

# 3. Refine Novelty Classification
results[, Final_Novelty := "Putative Novel"]

# A. Check against our comprehensive known list (RSID)
results[RSID %in% known_loci$RSID, Final_Novelty := "Known Locus"]

# B. Check against the fetched GWAS Catalog list
catalog_rsids <- fread("references/gwas_catalog_sle_rsids.tsv")$RSID
results[RSID %in% catalog_rsids, Final_Novelty := "Known Locus"]

# C. Distance-based check against the known loci coordinates (from step16)
known_coords <- data.table(
  CHR = c(1, 1, 1, 1, 1, 2, 2, 2, 4, 5, 5, 5, 6, 6, 6, 7, 8, 8, 11, 16, 16, 19, 2), # added 2 for STAT4
  POS = c(161509020, 169800269, 184982348, 2446426, 2696414, 203065200, 111894798, 107518970,
          12288952, 10247038, 96395099, 151078585, 32760884, 31234592, 137922602, 128945562,
          10958061, 11493510, 35080191, 85933077, 11097715, 10349293, 191094763),
  Gene = c("PTPN22", "FCGR2A", "TNFSF4", "ARID5B", "HP", "STAT4", "PFDN5", "XKR6",
           "BANK1", "BLK", "PTTG1", "TNIP1", "HLA-DRB1", "HLA-DQA1", "TNFAIP3", "IRF5",
           "XKR6", "BLK", "CD44", "IRF8", "CLEC16A", "TYK2", "STAT4")
)

for (i in which(results$Final_Novelty == "Putative Novel")) {
  chr <- results$CHR[i]
  pos <- results$BP[i]
  matches <- known_coords[CHR == chr & abs(POS - pos) < 500000]
  if (nrow(matches) > 0) {
    results$Final_Novelty[i] <- paste0("Known (near ", matches$Gene[1], ")")
  }
}

# 4. Add LAVA Replication Info and Meta
mapping <- fread("results/top_loci_lava_mapping.tsv")
results <- merge(results, mapping[, .(RSID, LAVA_LOC)], by = "RSID", all.x = TRUE)

lava <- fread("results/lava_sle_results.csv")
results <- merge(results, lava[, .(LAVA_LOC = LOC, LAVA_rho = rho, LAVA_P = p)], by = "LAVA_LOC", all.x = TRUE)

# Locus is High-Confidence (Replicated) if LAVA_P is Bonferroni significant & rho > 0
results[, Replicated := FALSE]
results[!is.na(LAVA_P) & LAVA_P < (0.05/nrow(results)) & LAVA_rho > 0, Replicated := TRUE]

# For legacy sensitivity scripts
rep_data <- fread("results/sensitivity/SA3_replication_power.tsv")
results <- merge(results, rep_data[, .(RSID, P_rep)], by = "RSID", all.x = TRUE)

# 5. Add Gene Annotation (placeholder or from previous table)
prev_table <- fread("results/top_loci_summary_table.tsv", fill=TRUE)
results <- merge(results, prev_table[, .(RSID, Gene_Annot = Gene, Region)], by = "RSID", all.x = TRUE)

# 6. Final Clean Table
final <- results[, .(
  RSID, CHR, BP, P_meta = P, BETA,
  Novelty = Final_Novelty,
  Replicated,
  LAVA_rho,
  LAVA_P,
  P_rep,
  Gene = Gene_Annot,
  Region
)]

# Save
fwrite(final, "results/final_locus_summary_table.tsv", sep = "\t")
# Overwrite the old top_loci table with refined data
fwrite(final, "results/top_loci_summary_table.tsv", sep = "\t")

cat("\n=== Final Summary Table Created ===\n")
cat("Total Loci:", nrow(final), "\n")
cat("Novel Loci:", sum(final$Novelty == "Putative Novel"), "\n")
cat("Known Loci:", sum(final$Novelty != "Putative Novel"), "\n")
cat("\nTop 5 Novel & Replicated Loci:\n")
print(head(final[Novelty == "Putative Novel" & Replicated == TRUE][order(P_meta)]))

cat("\nFinished.\n")
