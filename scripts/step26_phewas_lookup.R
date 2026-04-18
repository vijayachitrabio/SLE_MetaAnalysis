#!/usr/bin/env Rscript
# scripts/step26_phewas_lookup.R
# Refined PheWAS audit via GWAS Catalog REST API v2
# Logic: Strongest association per trait, with detailed metadata.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
  library(stringr)
})

DIR_RESULTS <- "results"
FILE_MASTER  <- file.path(DIR_RESULTS, "master_results_table.tsv")
FILE_OUT     <- file.path(DIR_RESULTS, "phewas_summary_refined.tsv")

cat("=== Starting Enhanced PheWAS Audit (v2 API + Core Fixes) ===\n")

if (!file.exists(FILE_MASTER)) stop("Master results table not found.")

# 1. Load High-Confidence RSIDs
master <- fread(FILE_MASTER)
target_rsids <- master[Final_Assessment == "HIGH CONFIDENCE", unique(RSID)]
cat("Found", length(target_rsids), "high-confidence RSIDs to audit.\n")

# 2. Refined Query Function
query_gwas_refined <- function(rsid) {
  cat(sprintf("  Processing %s... ", rsid))
  
  url <- "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
  # Use rs_id as the primary query param for v2
  res <- GET(url, query = list(rs_id = rsid, size = 100), add_headers(Accept = "application/json"))
  
  if (status_code(res) != 200) {
    cat(sprintf("Status %d. Skipped.\n", status_code(res)))
    return(NULL)
  }
  
  data <- fromJSON(content(res, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  
  if (!"_embedded" %in% names(data) || length(data$`_embedded`$associations) == 0) {
    cat("None found.\n")
    return(NULL)
  }
  
  assoc_list <- data$`_embedded`$associations
  
  # Process each association record
  extracted <- lapply(assoc_list, function(a) {
      # Trait extraction (Reported and EFO)
      reported_trait <- if(!is.null(a$reported_trait)) paste(a$reported_trait, collapse = "; ") else "Unknown"
      efo_id         <- if(length(a$efo_traits) > 0) a$efo_traits[[1]]$efo_id else NA
      efo_trait      <- if(length(a$efo_traits) > 0) a$efo_traits[[1]]$efo_trait else reported_trait

      # Beta / OR / Risk Allele
      beta          <- if(!is.null(a$beta) && a$beta != "-") as.numeric(a$beta) else NA
      odds_ratio    <- if(!is.null(a$or_per_copy_num)) as.numeric(a$or_per_copy_num) else NA
      risk_allele   <- if(length(a$snp_effect_allele) > 0) a$snp_effect_allele[[1]] else NA
      
      data.table(
          RSID = rsid,
          Study_ID = a$accession_id,
          Reported_Trait = reported_trait,
          EFO_Trait = efo_trait,
          EFO_ID = efo_id,
          P_value = a$p_value,
          Beta = beta,
          Odds_Ratio = odds_ratio,
          Risk_Allele = risk_allele,
          Pubmed_ID = a$pubmed_id
      )
  })
  
  dt <- rbindlist(extracted, fill = TRUE)
  
  # CORE FIX: Keep only the strongest association per unique Trait per RSID
  dt_cleaned <- dt %>%
      group_by(RSID, EFO_Trait) %>%
      arrange(P_value) %>%
      slice(1) %>%
      ungroup() %>%
      as.data.table()
  
  cat(sprintf("Found %d unique traits.\n", nrow(dt_cleaned)))
  return(dt_cleaned)
}

# 3. Execution
all_refined_results <- list()
for (rsid in target_rsids) {
  res <- query_gwas_refined(rsid)
  if (!is.null(res)) all_refined_results[[rsid]] <- res
  Sys.sleep(0.5)
}

final_df <- rbindlist(all_refined_results)

if (nrow(final_df) == 0) {
  cat("Zero associations found for refined audit.\n")
} else {
  # 4. Final Categorization
  autoimmune_keywords <- c("Lupus", "SLE", "Arthritis", "Scleroderma", "Sjogren", "Thyroiditis", 
                          "Diabetes", "Crohn", "Colitis", "Celiac", "Multiple sclerosis", "Autoimmune", "Psoriasis")
  
  final_df[, Category := ifelse(grepl(paste(autoimmune_keywords, collapse="|"), EFO_Trait, ignore.case = TRUE), 
                                 "Immune-Mediated", "Other Trait")]
  
  # 5. Save
  fwrite(final_df, FILE_OUT, sep = "\t")
  cat("Refined PheWAS audit saved to:", FILE_OUT, "\n")
  
  # Console preview
  cat("\n=== Pleiotropy Summary (Strongest Sig per Trait) ===\n")
  print(final_df[, .(N_Loci = .N), by = Category])
}
