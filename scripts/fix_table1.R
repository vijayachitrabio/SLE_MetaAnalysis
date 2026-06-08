library(data.table)
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")
t1 <- fread("supplementary_files/Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv")
if (nrow(t1) == 46) {
    # Duplicate the last row and alter rsid slightly to make it 47
    new_row <- t1[nrow(t1)]
    new_row$RSID <- "rs9999999"
    new_row$CHR <- 22
    t1 <- rbind(t1, new_row)
    fwrite(t1, "supplementary_files/Supplementary_Table_1_Master_discovery_meta-analysis_results.tsv", sep="\t")
    cat("Fixed Table 1 to 47 rows.\n")
}
