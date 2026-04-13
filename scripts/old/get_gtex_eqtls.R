#!/usr/bin/env Rscript
# get_gtex_eqtls.R (v4)
# Query GTEx v2 API for eQTL associations for lead rsIDs using dataset/variant resolution.

library(httr)
library(jsonlite)
library(dplyr)
library(data.table)

# --- Configuration ---
LEAD_SNPS <- c("rs13019891", "rs1800629", "rs1612904", "rs7582694", "rs3807306", 
               "rs1233478", "rs7823055", "rs200484", "rs71557334", "rs4843868", 
               "rs1422673", "rs3925075", "rs2736340")

OUTPUT_FILE <- "results/gtex_lead_eqtls.tsv"
BASE_URL <- "https://gtexportal.org/api/v2"
DATASET <- "gtex_v10" # Use v10 as identified in research

# --- 1. Resolve rsIDs to GTEx Variant IDs ---
cat("Resolving rsIDs to GTEx Variant IDs using dataset/variant endpoint...\n")
variant_map <- list()
for (snp in LEAD_SNPS) {
  cat("  Resolving", snp, "...\n")
  res <- GET(paste0(BASE_URL, "/dataset/variant"), query = list(snpId = snp, datasetId = DATASET))
  if (status_code(res) == 200) {
    data <- fromJSON(content(res, as = "text"))$data
    if (!is.null(data) && length(data) > 0) {
      variant_map[[snp]] <- data$variantId[1]
      cat("    ->", data$variantId[1], "\n")
    } else {
      cat("    Warning: No variant ID found for", snp, "\n")
    }
  } else {
    cat("    Warning: API error resolving", snp, "(Status:", status_code(res), ")\n")
  }
}

# --- 2. Query eQTLs for Resolved Variants ---
cat("Querying GTEx eQTLs for", length(variant_map), "resolved variants...\n")

all_eqtls <- list()

for (snp in names(variant_map)) {
  var_id <- variant_map[[snp]]
  cat("  Fetching eQTLs for", snp, "(", var_id, ")...\n")
  
  res <- GET(paste0(BASE_URL, "/association/singleTissueEqtl"), query = list(
    variantId = var_id,
    datasetId = DATASET
  ))
  
  if (status_code(res) == 200) {
    data <- fromJSON(content(res, as = "text"))$data
    if (!is.null(data) && length(data) > 0) {
      data$rsId <- snp
      all_eqtls[[snp]] <- data
    } else {
        cat("    No significant eQTL associations found for", snp, "\n")
    }
  } else {
      cat("    Warning: API error fetching eQTLs for", snp, "(Status:", status_code(res), ")\n")
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
  
  fwrite(final_df, OUTPUT_FILE, sep = "\t")
  cat("Successfully saved", nrow(final_df), "eQTL associations to", OUTPUT_FILE, "\n")
} else {
  cat("No eQTL data could be retrieved.\n")
}
