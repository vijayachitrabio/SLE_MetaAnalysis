library(openxlsx)
library(data.table)

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis/supplementary_files")

wb <- createWorkbook()
addWorksheet(wb, "Title Page")

# Define EXACT Titles (1 to 9 now!)
toc <- data.frame(
  "Table" = paste("Supplementary Table", 1:9),
  "Title" = c(
    "Genome-wide significant loci from the discovery meta-analysis",
    "Heterogeneity and cohort-level statistics for prioritized variants",
    "Consolidated LAVA prioritization results",
    "Bayesian colocalization results",
    "Spanish replication results for prioritized variants",
    "Pathway enrichment results",
    "PheWAS summary of prioritized SLE-associated variants",
    "Therapeutic annotation of prioritized SLE-associated genes",
    "Statistical fine-mapping (SuSiE) of prioritized SLE-associated loci"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

headerStyle <- createStyle(fontSize = 12, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD", border="TopBottom", borderColour = "#4F81BD")
writeDataTable(wb, "Title Page", x = toc, colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleMedium2", headerStyle = headerStyle)
setColWidths(wb, "Title Page", cols = 1:2, widths = c(25, 70))

boldStyle <- createStyle(textDecoration = "bold")
addStyle(wb, "Title Page", style = boldStyle, rows = 2:10, cols = 1, gridExpand = TRUE)

# Added SuSiE and PIP to abbreviations
abbrev_data <- data.frame(
  "Abbreviation" = c("LAVA", "PP4", "FDR", "LDSC", "PheWAS", "eQTL", "SuSiE", "PIP"),
  "Definition" = c("Local Analysis of Variant Association", "Posterior Probability of hypothesis 4 (colocalization)", "False Discovery Rate", "Linkage Disequilibrium Score Regression", "Phenome-Wide Association Study", "expression Quantitative Trait Locus", "Sum of Single Effects", "Posterior Inclusion Probability"),
  stringsAsFactors = FALSE, check.names = FALSE
)
# shifted down due to extra row in TOC
writeDataTable(wb, "Title Page", x = abbrev_data, startRow = 13, colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleLight9")
setColWidths(wb, "Title Page", cols = 1:2, widths = c(20, 80))
writeData(wb, "Title Page", x = "Abbreviations", startRow = 12)
addStyle(wb, "Title Page", style = boldStyle, rows = 12, cols = 1)

files <- c(
  "Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv",
  "Supplementary_Table_2_Quality-control_metrics_for_discovery_meta-analysis.tsv",
  "Supplementary_Table_3_Consolidated_LAVA_prioritization_results.tsv",
  "Supplementary_Table_4_Bayesian_colocalization_results.tsv",
  "Supplementary_Table_5_Spanish_replication_results_for_prioritized_variants.tsv",
  "Supplementary_Table_6_Pathway_enrichment_results.tsv",
  "Supplementary_Table_7_Cross-trait_genetic_correlation_PheWAS_summary.tsv",
  "Supplementary_Table_8_Drug_and_therapeutic_target_annotation_summary.tsv",
  "Supplementary_Table_9_SuSiE_statistical_fine-mapping_results.tsv"
)

titleStyle <- createStyle(fontSize = 14, textDecoration = "bold")

for (i in 1:9) {
  sheet_name <- paste("Supp_Table", i, sep="_")
  addWorksheet(wb, sheet_name)
  
  if (file.exists(files[i])) {
    dt <- fread(files[i])
    full_title <- paste0("Supplementary Table ", i, ". ", toc$Title[i])
    writeData(wb, sheet_name, x = full_title, startRow = 1, startCol = 1)
    addStyle(wb, sheet_name, style = titleStyle, rows = 1, cols = 1)
    writeDataTable(wb, sheet_name, x = dt, startRow = 3, tableStyle = "TableStyleLight9", withFilter = FALSE)
    setColWidths(wb, sheet_name, cols = 1:ncol(dt), widths = "auto")
  } else {
    writeData(wb, sheet_name, "File not found.", startRow = 3)
  }
}

saveWorkbook(wb, "Master_Supplementary_Tables.xlsx", overwrite = TRUE)
cat("Excel file created successfully.\n")
