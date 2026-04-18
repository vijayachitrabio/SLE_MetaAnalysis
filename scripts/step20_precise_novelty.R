#!/usr/bin/env Rscript
# scripts/step20_precise_novelty.R
# Precise novelty identification using ieugwasr (LD clumping) and GWAS Catalog

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ieugwasr)
})

cat("=== Precise Locus & Novelty Analysis ===\n")

# 1. Load Meta-analysis results
cat("Loading meta-analysis results...\n")
discovery <- fread("results/discovery_meta_results.tsv")
sig_snps <- discovery[P_meta < 5e-8, .(rsid = RSID, pval = P_meta)]

cat("Found", nrow(sig_snps), "genome-wide significant SNPs.\n")

# 2. LD Clumping using ieugwasr (OpenGWAS LD reference)
cat("Performing LD clumping (r2=0.1, kb=500)...\n")
# ieugwasr::ld_clump requires a data frame with rsid and pval
clumped <- ld_clump(sig_snps, clump_r2 = 0.1, clump_kb = 500)

cat("Independent lead SNPs after LD clumping:", nrow(clumped), "\n")

# 3. Merge with full stats to get BETA, CHR, BP
lead_stats <- merge(clumped, 
                    discovery[, .(RSID, CHR, BP, EA, OA, BETA = BETA_meta)], 
                    by.x = "rsid", by.y = "RSID")

# 4. Load GWAS Catalog Known SNPs
cat("Loading GWAS Catalog SLE SNPs...\n")
known_rsids <- fread("references/gwas_catalog_sle_rsids.tsv")$RSID

# 5. Novelty Classification
# A. Direct RSID match
lead_stats[, Known_by_RSID := rsid %in% known_rsids]

# B. Distance-based check (if not matching by RSID)
# We also want to check if it's "close" to a known SNP.
# Since we don't have all Catalog coordinates yet, we'll use a local list 
# of well-known SLE loci as a secondary check, OR query OpenGWAS for LD.

# Let's use the local list from step19b as a baseline
known_loci_coords <- data.table(
  CHR = c(1, 1, 1, 1, 2, 2, 2, 2, 4, 5, 5, 5, 6, 6, 6, 7, 8, 8, 11, 16, 16, 19),
  POS = c(161509020, 169800269, 184982348, 173282717, 203065200, 113072292, 137922602, 111894798, 
           102423596, 10247038, 96395099, 151078585, 32760884, 31234592, 137922602, 128945562,
           10958061, 11493510, 35080191, 85933077, 11097715, 10349293),
  Gene = c("PTPN22", "FCGR2A", "TNFSF4", "TNFSF4", "STAT4", "IL1F10", "TNFAIP3", "XKR6",
            "SIAE", "BLK", "PTTG1", "TNIP1", "HLA-DRB1", "HLA-DQA1", "TNFAIP3", "IRF5",
            "XKR6", "BLK", "CD44", "IRF8", "CLEC16A", "TYK2")
)

lead_stats[, Novelty := "Novel"]
lead_stats[Known_by_RSID == TRUE, Novelty := "Known (RSID match)"]

# Distance check for the rest
for (i in which(lead_stats$Novelty == "Novel")) {
  chr <- lead_stats$CHR[i]
  bp <- lead_stats$BP[i]
  
  match <- known_loci_coords[CHR == chr & abs(POS - bp) < 500000]
  if (nrow(match) > 0) {
    lead_stats$Novelty[i] <- paste0("Known region (", match$Gene[1], ")")
  }
}

# 6. Replication Status
cat("Adding replication info...\n")
rep_data <- fread("results/sensitivity/SA3_replication_power.tsv")
lead_stats <- merge(lead_stats, rep_data[, .(RSID, P_rep, BETA_disco)], 
                    by.x = "rsid", by.y = "RSID", all.x = TRUE)

n_tests <- nrow(lead_stats)
bonf_thresh <- 0.05 / n_tests

lead_stats[, Replication := "Not replicated"]
lead_stats[!is.na(P_rep) & P_rep < 0.05 & sign(BETA) == sign(BETA_disco), Replication := "Nominal"]
lead_stats[!is.na(P_rep) & P_rep < bonf_thresh, Replication := "Bonferroni"]

# 7. Final Output
final_table <- lead_stats[, .(
  SNP = rsid, CHR, BP, P_meta = pval, BETA_disco = BETA,
  Novelty, Replication, P_rep
)]

fwrite(final_table, "results/precise_novelty_results.tsv", sep = "\t")

cat("\nSummary of Novelty:\n")
print(table(final_table$Novelty))

cat("\nSummary of Replication:\n")
print(table(final_table$Replication))

cat("\nTruly Novel & Replicated loci:\n")
print(final_table[Novelty == "Novel" & Replication != "Not replicated"])

cat("\nDone. Results saved to results/precise_novelty_results.tsv\n")
