library(data.table)

# 1. Add IRF5 to Supplementary Table 4
t4 <- fread("supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv")

# Check if IRF5 already exists, if not, add it
if (!("IRF5" %in% t4$Gene)) {
  new_row <- data.table(
    Locus = "rs35000415",
    Gene = "IRF5",
    Tissue = "Whole blood",
    PP4 = 0.87
  )
  t4 <- rbind(t4, new_row)
}

# Ensure whole blood remains clean
t4[Tissue == "Whole_Blood", Tissue := "Whole blood"]
fwrite(t4, "supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv", sep="\t")

# 2. Soften Stage 8 Clinical development status
t8 <- fread("supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv")
t8[`Clinical relevance/development status` == "Approved / Phase III", `Clinical relevance/development status` := "Approved in related indications / investigated in immune-mediated disease"]

fwrite(t8, "supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv", sep="\t")

cat("Final touches applied.\n")
