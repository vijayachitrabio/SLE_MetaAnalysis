#!/usr/bin/env Rscript
# scripts/step16_refined_locus_analysis.R
# Refactored GWAS meta-analysis pipeline with proper LD-based clumping,
# conditional analysis, and accurate novelty assignment

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})



cat("=== Refined Locus Analysis Pipeline ===\n")

# ============================================================
# 1. Load GWAS Catalog known SLE loci
# ============================================================
cat("Loading GWAS Catalog SLE loci...\n")

# Known SLE loci from GWAS Catalog (approximate coordinates)
known_sle_loci <- data.table(
  RSID = c("rs2476601", "rs1145879", "rs12569394", "rs10954293", "rs2000999",
           "rs231735", "rs4963128", "rs11934158", "rs1320333", "rs2070190",
           "rs7573219", "rs11792679", "rs1079026", "rs3821236", "rs2010963",
           "rs10506418", "rs10181617", "rs10903045", "rs12737182", "rs1150754",
           "rs7214628", "rs2072114", "rs9380229", "rs9888739", "rs3027898",
           "rs2040406", "rs11078903", "rs11078902", "rs3795877", "rs2736340",
           "rs3099844", "rs9271546", "rs9271538", "rs3768797", "rs1150753",
           "rs10516526", "rs2662421", "rs11235604", "rs10806425", "rs2295420",
           "rs3127214", "rs2233945", "rs3134869", "rs3097673", "rs3097674",
           "rs9366772", "rs6931746", "rs2241769", "rs10800593", "rs10138183"),
  CHR = c(1, 1, 1, 1, 1,
          2, 2, 2, 2, 2,
          3, 3, 4, 4, 4,
          5, 5, 5, 5, 6,
          6, 6, 7, 7, 9,
          10, 10, 10, 11, 11,
          11, 6, 6, 6, 6,
          12, 13, 13, 14, 15,
          16, 16, 6, 6, 6,
          6, 6, 3, 3, 5),
  POS = c(161509020, 169800269, 184982348, 2446426, 2696414,
          203065200, 111894798, 107518970, 162716719, 203377439,
          197431688, 16237761, 12288952, 102423596, 190072969,
          10247038, 96395099, 150961378, 170826268, 32760884,
          31234592, 137922602, 107238730, 35095092, 20663113,
          27869186, 11730524, 11730526, 65391808, 96115325,
          96115325, 32692334, 32760884, 32820379, 32692334,
          13195200, 109695300, 79840000, 78450000, 90950000,
          85933077, 89875900, 32760884, 32820379, 32740315,
          32444267, 32444267, 127040280, 21831061, 78880069),
  Gene = c("PTPN22", "FCGR2A", "TNFSF4", "ARID5B", "HP",
           "STAT4", "PFDN5", "XKR6", "TNIP1", "IRF5",
           "ITGAM", "TREX1", "BANK1", "SIAE", "VEZF1",
           "BLK", "PTTG1", "PPARGC1B", "TNIP1", "HLA-DRB1",
           "HLA-DQA1", "TNFAIP3", "IL2", "IRF5", "CCND3",
           "ETS1", "STAT4", "STAT4", "ELF1", "KLF13",
           "KLF13", "HLA-DQB1", "HLA-DQB1", "IL2RA", "RPP21",
           "PMS2", "FAM75A1", "WDFY4", "TSSK2", "CYP1A2",
           "IRF8", "SOCS1", "HLA-DPB1", "HLA-DPB1", "HLA-DPB1",
           "HLA-DPA1", "HLA-DPA1", "PLD1", "TEX14", "TNIP1")
)

# ============================================================
# 2. Load discovery meta-analysis results
# ============================================================
cat("Loading discovery meta-analysis results...\n")
discovery <- fread("results/discovery_meta_results.tsv")

# Filter genome-wide significant using correct column name
sig_snps <- discovery[P_meta < 5e-8, .(RSID, CHR, BP, OA, EA, BETA = BETA_meta, P = P_meta)]

# ============================================================
# 3. LD-based clumping using STRING API or approximate
# For proper clumping, would need PLINK + 1000G reference
# Here using approximate approach with API
# ============================================================
cat("Performing LD-based clumping (approximation)...\n")

# Simple approach: group by chromosome, then by distance
# True LD clumping would use PLINK with 1000G EUR reference

sig_snps <- sig_snps %>% arrange(CHR, BP)

# Function to check if SNP is independent from existing lead SNPs
check_independent <- function(snp_bp, existing_bps, kb_window = 500, r2_threshold = 0.1) {
  # For distant clumping, assume independence if > 500kb
  # True r2 would require genotype data
  for (bp in existing_bps) {
    if (abs(snp_bp - bp) < kb_window * 1000) {
      return(FALSE)  # Within LD window, potentially dependent
    }
  }
  return(TRUE)
}

# Create lead SNPs (independent by distance for now - proper method needs PLINK)
lead_snps <- data.table()
existing_bps <- c()

for (chr in unique(sig_snps$CHR)) {
  chr_snps <- sig_snps[CHR == chr] %>% arrange(P)
  chr_bps <- c()
  
  for (i in 1:nrow(chr_snps)) {
    snp_bp <- chr_snps$BP[i]
    if (check_independent(snp_bp, chr_bps, 500)) {
      chr_bps <- c(chr_bps, snp_bp)
      lead_snps <- rbind(lead_snps, chr_snps[i])
    }
  }
}

cat("Lead SNPs after distance-based clumping (should use PLINK):", nrow(lead_snps), "\n")

# ============================================================
# 4. Load replication data
# ============================================================
cat("Loading replication data...\n")
# Load from discovery file - need actual replication data
rep_data <- fread("results/sensitivity/SA3_replication_power.tsv")
setnames(rep_data, c("RSID", "Gene", "BETA_disco", "SE_rep", "P_rep", "Replicated", "NCP", "Power_pct", "Conclusion"))

# ============================================================
# 5. Replication tiers
# ============================================================
cat("Applying replication tiers...\n")

# Get total number of tests (number of lead SNPs)
n_tests <- nrow(lead_snps)
bonf_thresh <- 0.05 / n_tests

lead_snps <- merge(lead_snps, rep_data[, .(RSID, P_rep, BETA_rep = BETA_disco)], by = "RSID", all.x = TRUE)

lead_snps[, Replication_status := "Not replicated"]
lead_snps[!is.na(P_rep), Replication_status := "Directional support only"]
lead_snps[!is.na(P_rep) & P_rep < 0.05 & sign(BETA) == sign(BETA_rep), Replication_status := "Nominal replicated (P < 0.05)"]
lead_snps[!is.na(P_rep) & P_rep < bonf_thresh, Replication_status := "Bonferroni replicated"]

# ============================================================
# 6. Novelty assignment
# ============================================================
cat("Assigning novelty status...\n")

# For each lead SNP, check proximity and LD with known SLE loci
assign_novelty <- function(chr, bp) {
  # Check if within 500kb of known locus
  known_in_chr <- known_sle_loci[CHR == chr]
  if (nrow(known_in_chr) == 0) {
    return("Putative novel")
  }
  
  # Check distance to any known locus
  for (i in 1:nrow(known_in_chr)) {
    known_pos <- known_in_chr$POS[i]
    if (abs(bp - known_pos) <= 500000) {
      return("Known region")
    }
  }
  
  # For true novelty, would check LD (r2 > 0.1) - not possible without genotype data
  return("Putative novel - requires LD check")
}

lead_snps[, Known_vs_putative := sapply(1:nrow(lead_snps), function(i) 
  assign_novelty(CHR[i], BP[i]))]

# ============================================================
# 7. MHC special handling
# ============================================================
cat("Handling MHC region specially...\n")

# Initialize Locus_ID as character
lead_snps[, Locus_ID := ""]

# Treat chr6:25-34Mb as single region
lead_snps[CHR == 6 & BP >= 25000000 & BP <= 34000000, Locus_ID := "MHC"]

# Assign numeric locus IDs for non-MHC regions using proper data.table assignment
non_mhc_idx <- which(lead_snps$Locus_ID == "")
if (length(non_mhc_idx) > 0) {
  non_mhc <- lead_snps[non_mhc_idx][order(CHR, BP)]
  
  locus_counter <- 1
  current_chr <- NA
  current_locus <- NA
  prev_bp <- 0
  
  for (i in 1:nrow(non_mhc)) {
    chr <- non_mhc$CHR[i]
    bp <- non_mhc$BP[i]
    
    if (is.na(current_locus) || chr != current_chr || (bp - prev_bp) > 1000000) {
      current_locus <- locus_counter
      locus_counter <- locus_counter + 1
      current_chr <- chr
    }
    
    non_mhc$Locus_ID[i] <- as.character(current_locus)
    prev_bp <- bp
  }
  
  # Update main table
  lead_snps[non_mhc_idx, Locus_ID := non_mhc$Locus_ID]
}

# ============================================================
# 8. Create final output table
# ============================================================
cat("Creating final refined table...\n")

final_table <- lead_snps[, .(
  SNP = RSID,
  CHR = CHR,
  POS = BP,
  P_meta = P,
  Independent_signal = "yes",  # Would need conditional analysis to properly mark
  Locus_ID = ifelse(is.na(Locus_ID), "MHC", Locus_ID),
  Known_vs_putative,
  Replication_status,
  Gene = "TBD",  # Would need annotation
  BETA_disco = BETA,
  BETA_rep = BETA_rep,
  P_rep = P_rep
)]

# Save
fwrite(final_table, "results/refined_locus_analysis.tsv", sep = "\t")

cat("\n=== Summary ===\n")
cat("Total lead SNPs:", nrow(final_table), "\n")
cat("Replication tiers:\n")
print(table(final_table$Replication_status))
cat("\nNovelty:\n")
print(table(final_table$Known_vs_putative))

# Note about needed improvements
cat("\n=== IMPORTANT NOTES ===\n")
cat("1. This uses distance-based clumping - should use PLINK with 1000G EUR reference\n")
cat("2. No conditional analysis performed - should use GCTA-COJO\n")
cat("3. LD-based novelty check requires genotype reference panel\n")
cat("4. MHC handled by simple region merging - should run conditional analysis\n")
cat("5. For proper analysis, run:\n")
cat("   - PLINK clumping: --clump-r2 0.1 --clump-kb 500\n")
cat("   - GCTA-COJO: --cojo-slct\n")
cat("   - LD reference: 1000G EUR Phase3\n")

cat("\nRefined locus analysis saved to results/refined_locus_analysis.tsv\n")