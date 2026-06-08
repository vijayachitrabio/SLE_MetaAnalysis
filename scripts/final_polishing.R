library(data.table)

# 1. Strictly align CLIC1 row in Supplementary Table 8
t8 <- fread("supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv")

t8[`Gene/target` == "CLIC1", `Pathway` := "Ion channel / inflammatory regulation"]
t8[`Gene/target` == "CLIC1", `Genetic/regulatory support` := "Discovery + LAVA + colocalization"]
t8[`Gene/target` == "CLIC1", `Representative agents/pathway examples` := "None validated in SLE"]
t8[`Gene/target` == "CLIC1", `Clinical relevance/development status` := "Preclinical / functional validation required"]
t8[`Gene/target` == "CLIC1", `Evidence category` := "Emerging candidate"]
t8[`Gene/target` == "CLIC1", `Notes/caveats` := "ρLAVA = 0.621; PP4 = 0.93 in spleen and 0.88 in whole blood"]

# Save back Table 8
fwrite(t8, "supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv", sep="\t")

# 2. In Supplementary Table 4, ensure "Whole_Blood" is strictly changed to "Whole blood"
t4 <- fread("supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv")
t4[Tissue == "Whole_Blood", Tissue := "Whole blood"]
fwrite(t4, "supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv", sep="\t")

cat("Polishing completed.\n")
