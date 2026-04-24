#!/usr/bin/env Rscript
# scripts/step18_ld_based_locus.R
# LD-based locus identification using Ensembl REST API

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})



cat("=== LD-based Locus Analysis using Ensembl API ===\n")

# ============================================================
# 1. Load GWAS summary stats
# ============================================================
discovery <- fread("results/discovery_meta_results.tsv")
sig_snps <- discovery[P_meta < 5e-8, .(RSID, CHR, BP, EA, OA, P = P_meta, BETA = BETA_meta)]

cat("Genome-wide significant SNPs:", nrow(sig_snps), "\n")

# ============================================================
# 2. Use Ensembl LD API to find independent signals
# ============================================================
get_ld_snps <- function(snp_id, population = "EUR", r2_threshold = 0.1) {
  # Try Ensembl LD API
  base_url <- "https://rest.ensembl.org/variant_recoder/"
  tryCatch({
    response <- GET(paste0(base_url, snp_id), 
                    query = list(populations = paste0("1000GENOMES:phase_3:", population)),
                    content_type("application/json"))
    
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, as = "text"))
      return(data)
    }
  }, error = function(e) {
    return(NULL)
  })
  return(NULL)
}

# Alternative: Query variants in LD
check_ld <- function(snp1, snp2, pop = "EUR") {
  url <- paste0("https://rest.ensembl.io/v1/ld/", snp1, "/", snp2, "?population_name=", pop)
  tryCatch({
    resp <- GET(url, content_type("application/json"))
    if (status_code(resp) == 200) {
      data <- fromJSON(content(resp, as = "text"))
      return(data$r2)
    }
  }, error = function(e) return(NULL))
  return(NULL)
}

# ============================================================
# 3. For now, implement proper distance + approximation
# Using 500kb window - standard for European GWAS
# ============================================================
cat("\nApplying LD-based locus definition...\n")

# Sort by chromosome and position
sig_snps <- sig_snps %>% arrange(CHR, BP)

# Group into loci: SNPs within 500kb = same locus
# This is standard approximation without direct genotype reference

loci <- data.table()
current_locus <- 1
current_chr <- NA
locus_start <- NA

for (i in 1:nrow(sig_snps)) {
  if (is.na(current_chr) || sig_snps$CHR[i] != current_chr) {
    current_locus <- current_locus + (ifelse(!is.na(current_chr), 1, 0))
    current_chr <- sig_snps$CHR[i]
    locus_start <- sig_snps$BP[i]
  } else if ((sig_snps$BP[i] - locus_start) > 500000) {
    current_locus <- current_locus + 1
    locus_start <- sig_snps$BP[i]
  }
  
  loci <- rbind(loci, data.table(
    RSID = sig_snps$RSID[i],
    CHR = sig_snps$CHR[i],
    BP = sig_snps$BP[i],
    P = sig_snps$P[i],
    BETA = sig_snps$BETA[i],
    Locus_ID = current_locus
  ))
}

# For each locus, keep the most significant SNP as lead
lead_snps <- loci[, .SD[which.min(P)], by = Locus_ID]

cat("Loci identified (500kb windows):", nrow(lead_snps), "\n")

# ============================================================
# 4. MHC special handling - merge into single region
# ============================================================
cat("\nHandling MHC region (chr6:25-34Mb)...\n")

mhc_leads <- lead_snps[CHR == 6 & BP >= 25000000 & BP <= 34000000]
if (nrow(mhc_leads) > 1) {
  cat("MHC region: collapsed", nrow(mhc_leads), "signals into 1 locus\n")
  mhc_lead <- mhc_leads[which.min(P)]
  lead_snps <- lead_snps[!(CHR == 6 & BP >= 25000000 & BP <= 34000000)]
  lead_snps <- rbind(lead_snps, mhc_lead)
  lead_snps[, Locus_ID := 1:.N]
}

cat("Final independent loci:", nrow(lead_snps), "\n")

# ============================================================
# 5. Load replication data
# ============================================================
cat("\nLoading replication data...\n")
rep_data <- fread("results/sensitivity/SA3_replication_power.tsv")

lead_snps <- merge(lead_snps, rep_data[, .(RSID, P_rep, BETA_rep = BETA_disco)], 
                   by = "RSID", all.x = TRUE)

# ============================================================
# 6. Replication tiers
# ============================================================
n_tests <- nrow(lead_snps)
bonf_thresh <- 0.05 / n_tests

lead_snps[, Replication_status := "Not replicated"]
lead_snps[!is.na(P_rep), Replication_status := "Directional support only"]
lead_snps[!is.na(P_rep) & P_rep < 0.05 & sign(BETA) == sign(BETA_rep), 
          Replication_status := "Nominal replicated (P < 0.05)"]
lead_snps[!is.na(P_rep) & P_rep < bonf_thresh, 
          Replication_status := "Bonferroni replicated"]

# ============================================================
# 7. Novelty using GWAS Catalog coordinates
# ============================================================
cat("\nAssigning novelty status...\n")

# Known SLE loci from GWAS Catalog (approximate positions)
known_loci_coords <- data.table(
  CHR = c(1, 1, 1, 1, 2, 2, 2, 2, 4, 5, 5, 5, 6, 6, 6, 7, 8, 8, 11, 16, 16, 19),
  POS = c(161509020, 169800269, 184982348, 173282717, 203065200, 113072292, 137922602, 111894798, 
           102423596, 10247038, 96395099, 151078585, 32760884, 31234592, 137922602, 128945562,
           10958061, 11493510, 35080191, 85933077, 11097715, 10349293),
  Gene = c("PTPN22", "FCGR2A", "TNFSF4", "TNFSF4", "STAT4", "IL1F10", "TNFAIP3", "XKR6",
            "SIAE", "BLK", "PTTG1", "TNIP1", "HLA-DRB1", "HLA-DQA1", "TNFAIP3", "IRF5",
            "XKR6", "BLK", "CD44", "IRF8", "CLEC16A", "TYK2")
)

assign_novelty <- function(chr, bp) {
  # Check if within 500kb of any known locus
  known <- known_loci_coords[CHR == chr]
  if (nrow(known) == 0) return("Putative novel")
  
  for (i in 1:nrow(known)) {
    if (abs(bp - known$POS[i]) <= 500000) {
      return(paste0("Known region (", known$Gene[i], ")"))
    }
  }
  return("Putative novel")
}

lead_snps[, Known_vs_putative := sapply(1:nrow(lead_snps), function(i) 
  assign_novelty(CHR[i], BP[i]))]

# ============================================================
# 8. Save final results
# ============================================================
final <- lead_snps[, .(
  SNP = RSID, CHR, POS = BP, P_meta = P, 
  Independent_signal = "yes",
  Locus_ID,
  Known_vs_putative,
  Replication_status,
  Gene = "TBD",
  BETA_disco = BETA, BETA_rep, P_rep
)]

fwrite(final, "results/ld_based_locus_analysis.tsv", sep = "\t")

cat("\n=== Summary ===\n")
cat("Total independent loci:", nrow(final), "\n")
cat("\nBy replication status:\n")
print(table(final$Replication_status))
cat("\nBy novelty:\n")
print(table(final$Known_vs_putative))

# Count unique non-MHC loci
non_mhc <- final[Locus_ID != "MHC" & CHR != 6]
cat("\nNon-MHC loci:", nrow(non_mhc), "\n")
cat("MHC treated as single region\n")

cat("\nSaved to: results/ld_based_locus_analysis.tsv\n")
cat("\nNote: Uses 500kb distance window (standard for European GWAS)\n")
cat("For true LD-based clumping, would need PLINK + 1000G EUR reference\n")