#!/usr/bin/env Rscript
# scripts/step8_sensitivity.R
# Sensitivity analyses for SLE GWAS Meta-Analysis (North vs South European)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

RES <- "results"
SENS_DIR <- file.path(RES, "sensitivity")
FIG <- "figures"
dir.create(SENS_DIR, showWarnings = FALSE)

cat("=== Checking Sensitivity Analyses ===\n")

# 1. Load Data
# We need to perform DerSimonian-Laird Random Effects vs Fixed Effects for Top Loci.
# Since we only saved `BETA_meta`, `SE_meta`, we need `BETA_bentham`, `SE_bentham` etc.
meta_results <- fread(file.path(RES, "discovery_meta_results.tsv"))
top_loci <- fread(file.path(RES, "top_loci_summary_table.tsv"))

# Merge Full Co-hort data for Top 57 loci
top_full <- merge(top_loci, meta_results[, .(RSID, BETA_bentham, SE_bentham, BETA_finngen, SE_finngen)], by="RSID")

cat("--- SA-1: Random-effects (DerSimonian-Laird) vs Fixed-effects ---\n")
dl_random <- top_full[, {
  b  <- c(BETA_bentham, BETA_finngen)
  se <- c(SE_bentham, SE_finngen)
  valid <- !is.na(b) & !is.na(se)
  b  <- b[valid]
  se <- se[valid]
  w  <- 1 / se^2
  W  <- sum(w)
  
  # Fixed
  b_fixed <- sum(w * b) / W
  Q  <- sum(w * (b - b_fixed)^2)
  k  <- length(b)
  df <- k - 1
  
  c_dl <- W - sum(w^2) / W
  tau2 <- max(0, (Q - df) / c_dl)
  
  # Random
  w_re  <- 1 / (se^2 + tau2)
  b_re  <- sum(w_re * b) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))
  z_re  <- b_re / se_re
  p_re  <- 2 * pnorm(-abs(z_re))
  
  .(
    Gene = Gene,
    P_fixed = P_disco,
    I2 = I2,
    tau2 = tau2,
    BETA_random = b_re,
    P_random = p_re,
    Substantial_Diff = abs(exp(b_re) - exp(b_fixed)) / exp(b_fixed) > 0.15
  )
}, by = .(RSID, CHR, BP)]

dl_random[, Note := fcase(
  Substantial_Diff == TRUE & tau2 > 0, "Random and fixed effects diverge",
  tau2 == 0, "No between-study variance",
  default = "Consistent"
)]

fwrite(dl_random, file.path(SENS_DIR, "SA1_random_vs_fixed.tsv"), sep="\t")
cat(sprintf("Found %d loci where random differs substantially from fixed.\n", sum(dl_random$Substantial_Diff)))

cat("--- SA-2: HLA-region signal independence ---\n")
hla_snps <- top_loci[CHR == 6, ]
hla_snps <- hla_snps[order(BP)]

pairs <- expand.grid(i = seq_len(nrow(hla_snps)), j = seq_len(nrow(hla_snps))) %>% filter(i < j)
if (nrow(pairs) > 0) {
  hla_pairs <- pairs %>% mutate(
    SNP_A = hla_snps$RSID[i], SNP_B = hla_snps$RSID[j],
    Gene_A = hla_snps$Gene[i], Gene_B = hla_snps$Gene[j],
    BP_A = hla_snps$BP[i], BP_B = hla_snps$BP[j],
    Distance_kb = abs(BP_A - BP_B) / 1000
  ) %>% mutate(Independent = Distance_kb >= 500) # Our pruning window was 500kb
  
  fwrite(hla_pairs, file.path(SENS_DIR, "SA2_HLA_pairwise.tsv"), sep="\t")
  cat(sprintf("All HLA pairs are independent (>500kb): %s\n", all(hla_pairs$Independent)))
}

cat("--- SA-3: Replication Power in Southern Cohort ---\n")
n_testable <- 57
alpha_bonf <- 0.05 / n_testable
z_bonf <- qnorm(1 - alpha_bonf / 2)

# Power depends on the SE in the replication cohort.
# SE_rep ~ abs(BETA_rep / qnorm(P_rep / 2))
top_loci[, SE_rep := abs(BETA_rep / qnorm(P_rep / 2))]
rep_snps <- top_loci[!is.na(SE_rep)]

power_dt <- rep_snps[, .(
  RSID, Gene, BETA_disco, SE_rep, P_rep, Replicated
)]
power_dt[, NCP := abs(BETA_disco) / SE_rep]
power_dt[, Power_pct := (pnorm(NCP - z_bonf) + pnorm(-NCP - z_bonf)) * 100]

power_dt[, Conclusion := fcase(
  Replicated == TRUE, "Replicated",
  Power_pct < 50 & !Replicated, "Underpowered",
  Power_pct >= 50 & !Replicated, "Adequately powered - non-replicating"
)]

fwrite(power_dt, file.path(SENS_DIR, "SA3_replication_power.tsv"), sep="\t")
cat(sprintf("Adequately powered non-replicating variants: %d\n", sum(power_dt$Conclusion == "Adequately powered - non-replicating")))

cat("Sensitivity analysis complete.\n")
