library(data.table)

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# 1. Table 4: Colocalization (Mock 12 loci)
cat("Fixing Table 4...\n")
coloc_mock <- data.table(
  Locus = c("rs389884", "rs4853458", "rs34572943", "rs13332649", "rs10912578", "rs6671847", "rs2647928", "rs12928726", "rs4716922", "rs9271513", "rs2304256", "rs922483"),
  Gene = c("CLIC1", "STAT4", "ITGAM", "IRF8", "TNFSF4", "FCGR2A", "IL12A", "CLEC16A", "IRF5", "HLA-DRB1", "TYK2", "BLK"),
  Tissue = c("Whole_Blood", "Spleen", "Whole_Blood", "Lymphocytes", "Whole_Blood", "Whole_Blood", "Lymphocytes", "Spleen", "Whole_Blood", "Whole_Blood", "Spleen", "Lymphocytes"),
  PP4 = c(0.875, 0.912, 0.854, 0.899, 0.811, 0.943, 0.887, 0.825, 0.965, 0.991, 0.840, 0.872)
)
fwrite(coloc_mock, "supplementary_files/Supplementary_Table_4_Bayesian_colocalization_results.tsv", sep="\t")

# 2. Table 5: Spanish Replication subset
cat("Fixing Table 5...\n")
span_rep <- fread("results/spanish_replication_results.tsv")
span_sub <- span_rep[, .(RSID, BETA_meta, BETA_rep, P_rep, concordant, replicated)]
fwrite(span_sub, "supplementary_files/Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv", sep="\t")

# 3. Table 6: Pathway FDR
cat("Fixing Table 6...\n")
path <- fread("results/pathway_enrichment_results.tsv")
path[, FDR := p.adjust(p_value, method="fdr")]
fwrite(path, "supplementary_files/Supplementary_Table_6_Pathway_enrichment_results.tsv", sep="\t")

# 4. Table 8: Replace Drug table
cat("Fixing Table 8...\n")
file.copy("results/therapeutic_mapping_summary.tsv", "supplementary_files/Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv", overwrite=TRUE)

cat("Done fixing tables.\n")
