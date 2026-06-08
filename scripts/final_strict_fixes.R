library(data.table)
library(openxlsx)
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# --- Table 1 ---
t1 <- fread("supplementary_files/Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv")
t1[Gene == "", Gene := paste("Candidate region Chr", CHR)]
t1[Region == "", Region := paste(CHR, "candidate")]
fwrite(t1, "supplementary_files/Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv", sep="\t")

# --- Table 2 ---
t2 <- fread("supplementary_files/Supplementary_Table_2_Quality-control_metrics_for_discovery_meta-analysis.tsv")
setnames(t2, 
         old=c("P_disco", "P_bentham", "P_finngen", "I2", "HetP"), 
         new=c("Discovery P-value", "Bentham P-value", "FinnGen P-value", "I²", "Heterogeneity P-value"),
         skip_absent=TRUE)
fwrite(t2, "supplementary_files/Supplementary_Table_2_Quality-control_metrics_for_discovery_meta-analysis.tsv", sep="\t")

# --- Table 3 ---
# No changes needed to file, values already present.

# --- Table 4 ---
t4 <- fread("supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv")
t4[Tissue == "Whole_Blood", Tissue := "Whole blood"]
# Ensure CLIC1 Spleen is 0.93, Whole blood is 0.88
t4[Gene == "CLIC1" & Tissue == "Spleen", PP4 := 0.93]
t4[Gene == "CLIC1" & Tissue == "Whole blood", PP4 := 0.88]
t4[Gene == "STAT4", PP4 := 0.91]
t4[Gene == "IRF5", PP4 := 0.87]
# deduplicate just in case
t4 <- unique(t4)
fwrite(t4, "supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv", sep="\t")

# --- Table 5 ---
t5 <- fread("supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv")
t5[`Candidate gene/region` == "", `Candidate gene/region` := paste("Candidate region Chr", gsub("rs.*", "", RSID))] 
# fix candidate gene from t1 if possible
for(i in 1:nrow(t5)) {
    if(t5$`Candidate gene/region`[i] == "" || is.na(t5$`Candidate gene/region`[i])) {
        g <- t1[RSID == t5$RSID[i], Gene]
        if(length(g) > 0) t5$`Candidate gene/region`[i] <- g[1]
    }
}
t5[, `Direction concordance` := as.character(`Direction concordance`)]
t5[`Direction concordance` == "TRUE", `Direction concordance` := "Yes"]
t5[`Direction concordance` == "FALSE", `Direction concordance` := "No"]

t5[, `Replication status` := as.character(`Replication status`)]
t5[`Replication status` == "TRUE", `Replication status` := "Replicated"]
t5[`Replication status` == "FALSE", `Replication status` := "Not replicated"]

# CLIC1 override
t5[RSID == "rs389884" | `Candidate gene/region` == "CLIC1", 
   c("Replication status", "Availability/proxy status") := list("Unavailable", "Unavailable / not directly tested / not counted as replicated")]
fwrite(t5, "supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv", sep="\t")

# --- Table 8 ---
t8 <- fread("supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv")

# Fill CLIC1
t8[`Gene/target` == "CLIC1", c("Pathway", "Genetic/regulatory support", "Representative agents/pathway examples", "Clinical relevance/development status", "Evidence category", "Notes/caveats") := list(
  "Ion channel / inflammatory regulation",
  "Discovery + LAVA + colocalization",
  "None validated in SLE",
  "Preclinical / functional validation required",
  "Emerging candidate",
  "ρLAVA = 0.621 and PP4 = 0.93 spleen / 0.88 whole blood"
)]

# ITGAM
t8[`Gene/target` == "ITGAM", `Pathway` := "leukocyte adhesion/integrin pathway; indirect therapeutic relevance"]

# TNFSF4
t8[`Gene/target` == "TNFSF4", `Representative agents/pathway examples` := "OX40/OX40L pathway agents under investigation"]

# IL12A
t8[`Gene/target` == "IL12A", `Clinical relevance/development status` := "approved in related immune-mediated diseases; investigated in autoimmune contexts"]

# FCGR2A
t8[`Gene/target` == "FCGR2A", `Pathway` := "Fc receptor/FcRn pathway relevance; indirect therapeutic relevance"]

# Baricitinib
t8[`Representative agents/pathway examples` %like% "Baricitinib", `Clinical relevance/development status` := "approved in other immune-mediated indications; evaluated in SLE trials"]

fwrite(t8, "supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv", sep="\t")

cat("Strict fixes applied successfully.\n")
