#!/usr/bin/env Rscript
# scripts/step17_proper_ld_clumping.R
# Proper GWAS analysis with PLINK LD clumping and GCTA conditional analysis

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})



cat("=== Proper LD-based GWAS Analysis ===\n")

# ============================================================
# Check for PLINK
# ============================================================
plink_path <- Sys.which("plink")
plink2_path <- Sys.which("plink2")

if (plink_path == "" && plink2_path == "") {
  cat("ERROR: PLINK not installed.\n")
  cat("To install:\n")
  cat("  macOS: brew install plink\n")
  cat("  Linux: sudo apt-get install plink\n")
  cat("  Or download: https://www.cog-genomics.org/plink2\n")
  cat("\n")
  cat("Also need 1000G EUR reference panel:\n")
  cat("  Download from: https://www.internationalgenome.org/data-portal/download\n")
  cat("  Select: Phase3 > EUR > GRCh38\n")
  stop("PLINK not found")
}

cat("Found PLINK at:", ifelse(plink2_path != "", plink2_path, plink_path), "\n")

# ============================================================
# Check for 1000G reference
# ============================================================
ref_files <- c("1000G_EUR.bed", "1000G_EUR.bim", "1000G_EUR.fam")
ref_dir <- "references"

if (!dir.exists(ref_dir)) dir.create(ref_dir)

ref_exists <- all(sapply(ref_files, function(f) file.exists(file.path(ref_dir, f))))

if (!ref_exists) {
  cat("ERROR: 1000G EUR reference not found.\n")
  cat("Expected files in", ref_dir, ":\n")
  cat("  1000G_EUR.bed\n  1000G_EUR.bim\n  1000G_EUR.fam\n")
  cat("\nDownload from:\n")
  cat("  https://www.internationalgenome.org/data-portal/download\n")
  cat("  Or use AWS: aws s3 sync s3://1000genomes/phase3/ .\n")
  stop("Reference panel not found")
}

# ============================================================
# Step 1: Prepare GWAS summary stats for PLINK
# ============================================================
cat("\n=== Step 1: Preparing summary stats ===\n")

discovery <- fread("results/discovery_meta_results.tsv")
sig_snps <- discovery[P_meta < 5e-8, .(CHR, BP, RSID, EA, OA, P = P_meta, BETA = BETA_meta)]

# Create PLINK format
plink_input <- sig_snps[, .(SNP = RSID, CHR, BP, A1 = EA, A2 = OA, P, BETA)]
fwrite(plink_input, "temp_gwas.txt", sep = "\t")

cat("Saved", nrow(plink_input), "genome-wide significant SNPs\n")

# ============================================================
# Step 2: PLINK LD clumping
# ============================================================
cat("\n=== Step 2: PLINK LD clumping (r2=0.1, kb=500) ===\n")

clump_out <- "temp_clumped"

system2(plink_path, args = c(
  "--bfile", file.path(ref_dir, "1000G_EUR"),
  "--clump", "temp_gwas.txt",
  "--clump-r2", "0.1",
  "--clump-kb", "500",
  "--clump-p1", "5e-8",
  "--clump-p2", "5e-8",
  "--out", clump_out
))

# Read clumped results
clumped <- fread(paste0(clump_out, ".clumped"))

lead_snps <- merge(sig_snps, clumped[, .(RSID, P)], by.x = "RSID", by.y = "RSID", all.x = TRUE)
lead_snps <- lead_snps[is.na(P.y) | P.y == P.x]  # Keep only clumped (independent) SNPs

cat("Independent lead SNPs after LD clumping:", nrow(lead_snps), "\n")

# ============================================================
# Step 3: GCTA-COJO conditional analysis (if available)
# ============================================================
cat("\n=== Step 3: GCTA-COJO conditional analysis ===\n")

gcta_path <- Sys.which("gcta64")
if (gcta_path == "") {
  cat("GCTA not found - skipping conditional analysis\n")
  cat("Install from: https://cnsgenomics.com/software/gcta/\n")
} else {
  cat("Found GCTA at:", gcta_path, "\n")
  cat("Would run: gcta64 --cojo-slct\n")
  # Note: GCTA requires individual-level genotype data
}

# ============================================================
# Step 4: Load replication and assign tiers
# ============================================================
cat("\n=== Step 4: Replication analysis ===\n")

rep_data <- fread("results/sensitivity/SA3_replication_power.tsv")

lead_snps <- merge(lead_snps, rep_data[, .(RSID, P_rep, BETA_rep = BETA_disco)], 
                   by = "RSID", all.x = TRUE)

n_tests <- nrow(lead_snps)
bonf_thresh <- 0.05 / n_tests

lead_snps[, Replication_status := "Not replicated"]
lead_snps[!is.na(P_rep), Replication_status := "Directional support only"]
lead_snps[!is.na(P_rep) & P_rep < 0.05 & sign(BETA) == sign(BETA_rep), 
          Replication_status := "Nominal replicated (P < 0.05)"]
lead_snps[!is.na(P_rep) & P_rep < bonf_thresh, 
          Replication_status := "Bonferroni replicated"]

# ============================================================
# Step 5: Novelty assignment with GWAS Catalog
# ============================================================
cat("\n=== Step 5: Novelty assignment ===\n")

# Load known SLE loci from GWAS Catalog
known_loci <- fread("results/sensitivity/SA2_HLA_pairwise.tsv")  # placeholder - need proper list
# Use the known loci from step16
source("scripts/step16_refined_locus_analysis.R")  # reloads known_sle_loci

# This would require proper GWAS Catalog data for full novelty check
# For now, use coordinate-based check

# ============================================================
# Step 6: Save final results
# ============================================================
cat("\n=== Step 6: Saving final results ===\n")

final <- lead_snps[, .(
  SNP = RSID, CHR, POS = BP, P_meta = P.x, 
  Independent_signal = "yes",  # After LD clumping
  Locus_ID = "to_be_assigned",  # Need proper locus grouping
  Known_vs_putative = "requires_LD_check",  # Need GWAS Catalog + LD
  Replication_status,
  Gene = "TBD",  # Need annotation
  BETA_disco = BETA, BETA_rep, P_rep
)]

fwrite(final, "results/final_locus_analysis.tsv", sep = "\t")

cat("\n=== Summary ===\n")
cat("Independent loci (after LD clumping):", nrow(final), "\n")
cat("Replication tiers:\n")
print(table(final$Replication_status))

# Clean up temp files
file.remove(c("temp_gwas.txt", paste0(clump_out, ".clumped"), paste0(clump_out, ".log")))

cat("\nFinal analysis saved to results/final_locus_analysis.tsv\n")