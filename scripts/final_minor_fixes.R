library(data.table)
library(openxlsx)
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# --- Table 5 ---
t5 <- fread("supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv")
t5[RSID == "rs10912578", `Candidate gene/region` := "TNFSF4 / TRIM21-region"]
t5[RSID == "rs12928726", `Candidate gene/region` := "CLEC16A"]
t5[RSID == "rs13332649", `Candidate gene/region` := "IRF8-region"]
t5[RSID == "rs2647928", `Candidate gene/region` := "IL12A"]
t5[RSID == "rs34572943", `Candidate gene/region` := "ITGAM"]
t5[RSID == "rs4853458", `Candidate gene/region` := "STAT4"]
t5[RSID == "rs6671847", `Candidate gene/region` := "FCGR2A / STAT1-region"]
t5[RSID == "rs389884", `Candidate gene/region` := "CLIC1-region candidate signal"]
fwrite(t5, "supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv", sep="\t")

# --- Table 8 ---
t8 <- fread("supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv")
t8[`Gene/target` == "ITGAM", `Representative agents/pathway examples` := "Complement-pathway agents, e.g., eculizumab/C5, avacopan/C5aR; indirect relevance only."]
t8[`Gene/target` == "FCGR2A", `Representative agents/pathway examples` := "IVIG and FcRn-pathway agents; indirect Fc-receptor pathway relevance."]
fwrite(t8, "supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv", sep="\t")

# --- Table 1 and 3 Gene/Region consistency ---
t1 <- fread("supplementary_files/Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv")
t3 <- fread("supplementary_files/Supplementary_Table_3_Consolidated_LAVA_prioritization_results.tsv")

# Fix rs17849501 and rs41272536 which are in LAVA block 138 (HIGH CONFIDENCE)
# Based on final_lava_consolidated_loci.tsv, they had empty genes.
t1[RSID %in% c("rs17849501", "rs41272536") & Gene %like% "Candidate", Gene := "1q25 locus"]
t1[RSID %in% c("rs17849501", "rs41272536") & Region %like% "candidate", Region := "1q25"]
t3[RSID %in% c("rs17849501", "rs41272536") & Gene == "", Gene := "1q25 locus"]
t3[RSID %in% c("rs17849501", "rs41272536") & Region == "", Region := "1q25"]

fwrite(t1, "supplementary_files/Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv", sep="\t")
fwrite(t3, "supplementary_files/Supplementary_Table_3_Consolidated_LAVA_prioritization_results.tsv", sep="\t")

cat("Minor fixes applied successfully.\n")
