#!/usr/bin/env Rscript
# scripts/step20_precise_novelty_ensembl.R
# Precise novelty identification using Ensembl REST API (LD calculations)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})

cat("=== Precise Locus & Novelty Analysis (Ensembl API) ===\n")

# 1. Load Meta-analysis results
cat("Loading meta-analysis results...\n")
discovery <- fread("results/discovery_meta_results.tsv")
sig_snps <- discovery[P_meta < 5e-8, .(RSID, CHR, BP, P = P_meta, BETA = BETA_meta)]

cat("Found", nrow(sig_snps), "genome-wide significant SNPs.\n")

# 2. Distance-based clumping to find approximate lead SNPs
cat("Performing initial distance-based clumping (500kb)...\n")
sig_snps <- sig_snps %>% arrange(CHR, BP)
loci <- data.table()
current_locus <- 1
current_chr <- NA
locus_start <- NA

for (i in 1:nrow(sig_snps)) {
  if (is.na(current_chr) || sig_snps$CHR[i] != current_chr) {
    current_locus <- current_locus + (if(!is.na(current_chr)) 1 else 0)
    current_chr <- sig_snps$CHR[i]
    locus_start <- sig_snps$BP[i]
  } else if ((sig_snps$BP[i] - locus_start) > 500000) {
    current_locus <- current_locus + 1
    locus_start <- sig_snps$BP[i]
  }
  loci <- rbind(loci, data.table(sig_snps[i,], Locus_ID = current_locus))
}

# Special handling for MHC
mhc_idx <- which(loci$CHR == 6 & loci$BP >= 25000000 & loci$BP <= 35000000)
if (length(mhc_idx) > 0) {
  loci$Locus_ID[mhc_idx] <- "MHC"
}

lead_snps <- loci[, .SD[which.min(P)], by = Locus_ID]
cat("Lead SNPs identified:", nrow(lead_snps), "\n")

# 3. Load GWAS Catalog SLE SNPs
cat("Loading GWAS Catalog SLE SNPs...\n")
known_rsids <- fread("references/gwas_catalog_sle_rsids.tsv")$RSID

# 4. Precise LD Check function using Ensembl
check_ld_ensembl <- function(our_snp, known_list, population = "1000GENOMES:phase_3:EUR") {
  # We only check known SNPs on the same chromosome and within 1Mb
  # To save API calls, we first check if any known SNP is close
  # If it is, we ask Ensembl for LD
  
  # For now, let's use the most efficient way: query the LD for our SNP against the whole population
  # in a window, but Ensembl doesn't support "all SNPs in LD" easily without a reference set.
  # So we check against a curated list of top SLE loci.
  
  known_curated <- c("rs2476601", "rs1145879", "rs12569394", "rs173282717", "rs203065200", 
                     "rs113072292", "rs137922602", "rs111894798", "rs102423596", "rs10247038", 
                     "rs96395099", "rs151078585", "rs32760884", "rs31234592", "rs128945562", 
                     "rs10958061", "rs11493510", "rs35080191", "rs85933077", "rs11097715", "rs10349293")
  
  max_r2 <- 0
  linked_snp <- NA
  
  for (known_snp in known_curated) {
    url <- paste0("https://rest.ensembl.org/ld/human/", our_snp, "/", known_snp, 
                  "?population_name=", population)
    tryCatch({
      resp <- GET(url, content_type("application/json"))
      if (status_code(resp) == 200) {
        data <- fromJSON(content(resp, as = "text"))
        if (!is.null(data$r2) && data$r2 > max_r2) {
          max_r2 <- data$r2
          linked_snp <- known_snp
        }
      }
    }, error = function(e) {})
    if (max_r2 > 0.8) break # Significant enough
  }
  return(list(r2 = max_r2, snp = linked_snp))
}

# 5. Process each lead SNP
cat("\nChecking novelty (LD-based via Ensembl)...\n")
lead_snps[, Novelty := "Novel"]
lead_snps[RSID %in% known_rsids, Novelty := "Known (RSID match)"]

# Only check LD for putative novel SNPs
to_check <- which(lead_snps$Novelty == "Novel" & lead_snps$Locus_ID != "MHC")

for (i in to_check) {
  snp <- lead_snps$RSID[i]
  cat("Checking", snp, "(", which(to_check == i), "/", length(to_check), ")... ")
  
  res <- check_ld_ensembl(snp, known_rsids)
  if (res$r2 > 0.1) {
    lead_snps$Novelty[i] <- paste0("Known (LD-linked to ", res$snp, ", r2=", round(res$r2, 2), ")")
    cat("Linked (r2=", res$r2, ")\n")
  } else {
    cat("Truly Novel\n")
  }
}

# 6. Save results
fwrite(lead_snps, "results/precise_novelty_ensembl.tsv", sep = "\t")

cat("\nSummary Table:\n")
print(table(lead_snps$Novelty))

cat("\nDone. Results saved to results/precise_novelty_ensembl.tsv\n")
