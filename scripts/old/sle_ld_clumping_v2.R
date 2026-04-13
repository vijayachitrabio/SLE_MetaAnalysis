library(data.table)
library(dplyr)

RESULTS_DIR <- "../results"
dir.create(RESULTS_DIR, showWarnings = FALSE)

cat("=== Conservative Physical Distance Pruning v2 ===\n")
cat("NOTE: This script uses physical distance pruning (not LD-based clumping).\n")
cat("      Results are described as 'lead variants after conservative physical\n")
cat("      distance pruning' or 'putatively independent lead variants'.\n")
cat("      For true LD-based independence use sle_ld_clumping_plink.R.\n\n")

# Load meta-analysis results
meta_dt <- fread(file.path(RESULTS_DIR, "sle_meta_analysis_results.tsv"))
cat("Total SNPs in meta-analysis:", nrow(meta_dt), "\n")

# Filter for genome-wide significance (5e-8)
sig_snps <- meta_dt[P < 5e-8]
cat("Significant SNPs (P < 5e-8):", nrow(sig_snps), "\n")

if (nrow(sig_snps) == 0) {
    stop("No significant SNPs to clump.")
}

# Define clumping window sizes
# MHC/HLA region (Chr 6: 25-34 Mb) prone to long-range LD
WINDOW_HLA <- 1000000  # 1 MB
# Standard genome-wide LD window for Europeans
WINDOW_DEFAULT <- 500000  # 500 KB

# Sort by P-value (best signal first)
sig_snps <- sig_snps[order(P)]

lead_snps <- data.table()
remaining_snps <- copy(sig_snps)

while (nrow(remaining_snps) > 0) {
    # Pick the top SNP
    top_snp <- remaining_snps[1, ]
    lead_snps <- rbind(lead_snps, top_snp)
    
    # Define window for this SNP
    current_window <- ifelse(top_snp$CHR == 6 & top_snp$BP >= 25000000 & top_snp$BP <= 34000000, 
                             WINDOW_HLA, 
                             WINDOW_DEFAULT)
    
    # Remove SNPs within window on the same chromosome
    remaining_snps <- remaining_snps[!(CHR == top_snp$CHR & 
                                       BP >= top_snp$BP - current_window & 
                                       BP <= top_snp$BP + current_window)]
    
    cat(sprintf("Picked %s (Chr%s:%s, P=%.2e), excluded SNPs within %d kb.\n", 
                top_snp$SNP, top_snp$CHR, top_snp$BP, top_snp$P, current_window/1000))
}

cat("Total putatively independent lead variants (conservative physical pruning):", nrow(lead_snps), "\n")

# Save results
fwrite(lead_snps, file.path(RESULTS_DIR, "sle_meta_lead_snps_v2.tsv"), sep="\t")
# Overwrite main lead SNPs file for downstream scripts
fwrite(lead_snps, file.path(RESULTS_DIR, "sle_meta_lead_snps.tsv"), sep="\t")

cat("Lead SNPs saved to results/sle_meta_lead_snps.tsv\n")
print(lead_snps[, .(SNP, CHR, BP, BETA, SE, P)])

cat("\n=== Clumping complete ===\n")
