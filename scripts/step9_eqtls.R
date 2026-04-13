#!/usr/bin/env Rscript
# scripts/step9_eqtls.R
# Query GTEx v8 API for cis-eQTL evidence on replicated loci.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})

RES <- "results"
cat("=== Querying GTEx for eQTL Evidence ===\n")

# Load replicated loci
top_loci <- fread(file.path(RES, "top_loci_summary_table.tsv"))
rep_loci <- top_loci[Replicated == TRUE, ]

# Query GTEx v2 API for eQTL associations for replicated lead rsIDs

BASE_URL <- "https://gtexportal.org/api/v2"
DATASET <- "gtex_v10"

# --- 1. Resolve rsIDs to GTEx Variant IDs ---
cat("Resolving rsIDs to GTEx Variant IDs...\n")
variant_map <- list()
for (snp in rep_loci$RSID) {
  cat("  Resolving", snp, "...\n")
  res <- tryCatch(GET(paste0(BASE_URL, "/dataset/variant"), query = list(snpId = snp, datasetId = DATASET)), error=function(e) NULL)
  if (!is.null(res) && status_code(res) == 200) {
    data <- fromJSON(content(res, as = "text"))$data
    if (!is.null(data) && length(data) > 0) {
      variant_map[[snp]] <- data$variantId[1]
    } 
  }
}

# --- 2. Query eQTLs for Resolved Variants ---
cat("Querying GTEx eQTLs for", length(variant_map), "resolved variants in immune tissues...\n")
immune_tissues <- c("Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes")
all_eqtls <- list()

for (snp in names(variant_map)) {
  var_id <- variant_map[[snp]]
  
  res <- tryCatch(GET(paste0(BASE_URL, "/association/singleTissueEqtl"), query = list(
    variantId = var_id,
    datasetId = DATASET
  )), error=function(e) NULL)
  
  if (!is.null(res) && status_code(res) == 200) {
    data <- fromJSON(content(res, as = "text"))$data
    if (!is.null(data) && length(data) > 0) {
      # Filter to immune tissues and significance
      sig_eqtls <- data %>% 
        filter(tissueSiteDetailId %in% immune_tissues, pValue < 1e-5) %>%
        arrange(pValue) %>% head(5)
        
      if(nrow(sig_eqtls) > 0) {
        sig_eqtls$rsId <- snp
        all_eqtls[[snp]] <- sig_eqtls
      }
    }
  }
  Sys.sleep(0.3) 
}

# --- Combine and Save ---
if (length(all_eqtls) > 0) {
  final_df <- bind_rows(all_eqtls)
  # Basic data cleaning
  final_df <- final_df %>%
    rename(tissue = tissueSiteDetailId, pval = pValue, nes = nes) %>%
    select(rsId, variantId, geneSymbol, gencodeId, tissue, pval, nes)
  
  fwrite(final_df, file.path(RES, "eqtl_summary.tsv"), sep = "\t")
  cat("Successfully saved", nrow(final_df), "eQTL associations to results/eqtl_summary.tsv\n")
} else {
  cat("No significant immune eQTL data could be retrieved for the variants.\n")
}

