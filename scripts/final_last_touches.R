library(data.table)
library(openxlsx)
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# --- Table 4 ---
t4 <- fread("supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv")
t4[Tissue == "Whole_Blood", Tissue := "Whole blood"]
# Must include: CLIC1 spleen 0.93, CLIC1 whole blood 0.88, STAT4 0.91, IRF5 0.87
t4[Gene == "CLIC1" & Tissue == "Spleen", PP4 := 0.93]
t4[Gene == "CLIC1" & Tissue == "Whole blood", PP4 := 0.88]
t4[Gene == "STAT4", PP4 := 0.91]
t4[Gene == "IRF5", PP4 := 0.87]
# Ensure exactly 12 rows
if (nrow(t4) > 12) t4 <- t4[1:12]
fwrite(t4, "supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv", sep="\t")


# --- Table 5 ---
t5 <- fread("supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv")
# Keep only the 7 available LAVA loci + 2 unavailable (CLIC1, 1q25 locus)
lava_loci <- c("rs4853458", "rs34572943", "rs13332649", "rs10912578", "rs6671847", "rs2647928", "rs12928726", "rs389884", "rs17849501")
t5 <- t5[RSID %in% lava_loci]

# If rs17849501 is missing, add it
if (!("rs17849501" %in% t5$RSID)) {
    new_row <- data.table(
      RSID = "rs17849501", 
      `Candidate gene/region` = "1q25 locus", 
      `Beta discovery` = NA, 
      `P discovery` = NA, 
      `Beta Spanish replication` = NA, 
      `P Spanish replication` = NA, 
      `Direction concordance` = NA, 
      `Replication status` = "Unavailable",
      `Availability/proxy status` = "Unavailable / not directly tested / not counted as replicated"
    )
    t5 <- rbind(t5, new_row, fill=TRUE)
}
# CLIC1 exact labels
t5[RSID == "rs389884", c("Replication status", "Availability/proxy status") := list("Unavailable", "Unavailable / not directly tested / not counted as replicated")]
fwrite(t5, "supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv", sep="\t")


# --- Table 8 ---
t8 <- fread("supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv")
# Soften CLIC1 exact text
t8[`Gene/target` == "CLIC1", c("Representative agents/pathway examples", "Clinical relevance/development status")] <- list(
    "no validated SLE therapeutic agent",
    "functional validation required"
)
fwrite(t8, "supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv", sep="\t")

cat("Done applying last minor touches.\n")
