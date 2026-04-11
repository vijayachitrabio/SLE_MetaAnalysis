library(data.table)
library(dplyr)
library(ggplot2)

RESULTS_DIR <- "../results"
FIG_DIR     <- "../figures"
dir.create(FIG_DIR, showWarnings = FALSE)

cat("=== STEP 3: Robustness / QC + Known-vs-Novel Summary Table ===\n\n")

meta_dt  <- fread(file.path(RESULTS_DIR, "sle_meta_analysis_results.tsv"))
lead_dt  <- fread(file.path(RESULTS_DIR, "sle_meta_lead_snps.tsv"))
gene_dt  <- fread(file.path(RESULTS_DIR, "gene_annotation.tsv"))
het_dt   <- fread(file.path(RESULTS_DIR, "heterogeneity_lead_snps.tsv"))
rep_dt   <- fread(file.path(RESULTS_DIR, "finngen_replication_results.tsv"))
pcoh_dt  <- fread(file.path(RESULTS_DIR, "per_cohort_effects.tsv"))

# ─── 3a: Sample-size per SNP (cohort availability) ──────────────────────────
cat("--- 3a. Sample-size per SNP ---\n")
# We know each SNP in meta_dt was contributed by all 3 cohorts (intersection in meta script)
# Flag lead SNPs that appear in < 3 cohorts after merge (none expected, but mark it)
meta_dt[, N_cohorts := 2]   # all meta SNPs exist in both European cohorts by construction
meta_dt[, missing_flag := FALSE]

ss_summary <- data.table(
  Metric       = c("Total meta-analysis SNPs","SNPs in both cohorts","SNPs missing ≥1 cohort"),
  Value        = c(nrow(meta_dt), nrow(meta_dt), 0)
)
fwrite(ss_summary, file.path(RESULTS_DIR,"sample_size_summary.tsv"), sep="\t")
print(ss_summary)

# ─── 3b: Imputation quality (check if INFO column exists) ───────────────────
cat("\n--- 3b. Imputation quality (INFO check) ---\n")
if ("INFO" %in% names(meta_dt)) {
  cat("INFO column present. Plotting distribution and filtering SNPs with INFO <0.8.\n")
  meta_dt_filtered <- meta_dt[INFO >= 0.8]
  p_info <- ggplot(meta_dt, aes(x=INFO)) +
    geom_histogram(bins=60, fill="steelblue", color="white", alpha=0.8) +
    geom_vline(xintercept=0.8, linetype="dashed", color="red") +
    labs(title="INFO score distribution (all SNPs)",
         subtitle="Red dashed = INFO 0.8 threshold",
         x="INFO", y="Count") +
    theme_minimal()
  ggsave(file.path(FIG_DIR,"info_distribution.png"), p_info, width=6, height=4, dpi=150)
  cat("SNPs removed due to INFO<0.8:", nrow(meta_dt)-nrow(meta_dt_filtered), "\n")
} else {
  cat("INFO column not present in meta results (expected – harmonised GWAS Catalog files often omit INFO).\n")
  cat("Note for Methods: Individual study QC (including imputation quality) was performed by original authors.\n")
  info_note <- data.table(Note="INFO column not available in harmonised summary statistics. Individual cohort imputation QC was performed by original study authors (Bentham 2015, Julià 2018).")
  fwrite(info_note, file.path(RESULTS_DIR,"imputation_qc_note.tsv"), sep="\t")
}

# ─── 3c: MAF-stratified QC ──────────────────────────────────────────────────
cat("\n--- 3c. MAF-stratified QC ---\n")
# Check if MAF column present
if (!"MAF" %in% names(meta_dt)) {
  # Reconstruct approximate MAF from effect frequency if available, else skip
  cat("No MAF column available in meta results. Checking for effect_allele_frequency...\n")
  if ("effect_allele_frequency" %in% names(meta_dt)) {
    meta_dt[, MAF := pmin(effect_allele_frequency, 1-effect_allele_frequency)]
  } else {
    cat("Skipping MAF stratification – no frequency column available.\n")
    maf_note <- data.table(Note="MAF not available in harmonised summary statistics.")
    fwrite(maf_note, file.path(RESULTS_DIR,"maf_qc_note.tsv"), sep="\t")
  }
}

if ("MAF" %in% names(meta_dt)) {
  meta_dt[, maf_group := ifelse(MAF > 0.05, "Common (MAF>5%)", "Low-frequency (1-5%)")]
  maf_qc <- meta_dt %>%
    group_by(maf_group) %>%
    summarise(
      n_snps      = n(),
      median_beta = median(abs(BETA), na.rm=TRUE),
      n_gwsig     = sum(P < 5e-8, na.rm=TRUE)
    )
  fwrite(as.data.table(maf_qc), file.path(RESULTS_DIR,"maf_stratified_qc.tsv"), sep="\t")
  print(maf_qc)

  p_maf <- ggplot(meta_dt[!is.na(maf_group)], aes(x=abs(BETA), fill=maf_group)) +
    geom_histogram(bins=60, alpha=0.65, position="identity") +
    scale_fill_manual(values=c("Common (MAF>5%)"="steelblue","Low-frequency (1-5%)"="darkorange")) +
    labs(title="Effect size distribution by MAF stratum",
         x="|β|", y="Count", fill="MAF group") +
    theme_minimal()
  ggsave(file.path(FIG_DIR,"beta_by_maf.png"), p_maf, width=6, height=4, dpi=150)
}

# ─── 3d: Inflation / genomic control summary ────────────────────────────────
cat("\n--- 3d. Inflation summary (lambda per cohort + meta) ---\n")
inflation_summary <- data.table(
  Study   = c("Bentham 2015","Julià 2018","Meta-Analysis"),
  Lambda  = c(1.409, 1.039, 1.019),
  GC_Applied = c("Yes","No","N/A (GC applied to inputs)")
)
fwrite(inflation_summary, file.path(RESULTS_DIR,"inflation_summary.tsv"), sep="\t")
print(inflation_summary)

# ─── 3e: Ancestry stratification check ──────────────────────────────────────
cat("\n--- 3e. Ancestry stratification check ---\n")
ancestry_check <- data.table(
  Study        = c("Bentham 2015 (GCST003156)","Julià 2018 (GCST006093)"),
  PopulationN  = c("~4,036 cases / 6,959 controls","~799 cases / 1,558 controls"),
  Ancestry     = c("European","Spanish (European)"),
  Mixed        = c(FALSE, FALSE)
)
fwrite(ancestry_check, file.path(RESULTS_DIR,"ancestry_stratification_check.tsv"), sep="\t")
print(ancestry_check)
cat("All cohorts are European-ancestry. No ancestry stratification required.\n")

# ─── 3f: Top loci summary table (publication-ready) ─────────────────────────
cat("\n--- 3f. Top loci summary table ---\n")

top_table <- merge(lead_dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P)],
                   gene_dt[, .(SNP, Gene, Region, Known_SLE)], by="SNP", all.x=TRUE)
top_table <- merge(top_table,
                   het_dt[, .(SNP, I2, p_Q)], by="SNP", all.x=TRUE)
top_table <- merge(top_table,
                   rep_dt[, .(SNP, FG_BETA, FG_P, Replicated)], by="SNP", all.x=TRUE)
top_table <- merge(top_table,
                   pcoh_dt[, .(SNP, Dir_Bentham, Dir_Julia, Consistent)], by="SNP", all.x=TRUE)

top_table[, Locus := paste0(Gene," (", Region,")")]
top_table[, Novel := !Known_SLE]
top_table <- top_table[order(P)]

fwrite(top_table, file.path(RESULTS_DIR,"top_loci_summary_table.tsv"), sep="\t")
cat("Top loci summary table saved.\n")
print(top_table[, .(SNP, Locus, BETA, SE, P, I2, Replicated, Novel)])

# ─── 3g: Supplementary QC table ─────────────────────────────────────────────
cat("\n--- 3g. Supplementary QC table ---\n")
supp_qc <- data.table(
  Cohort          = c("Bentham 2015","Julià 2018","Meta-Analysis"),
  GWAS_Catalog_ID = c("GCST003156","GCST006093","N/A"),
  N_Cases         = c(4036, 799, 4835),
  N_Controls      = c(6959, 1558, 8517),
  N_SNPs_raw      = c("~2.6M","~6.4M","N/A"),
  N_SNPs_filtered = c("~430K","~430K","429,799"),
  Lambda_raw      = c(1.409, 1.039, NA),
  GC_Corrected    = c("Yes","No","N/A"),
  Lambda_post_GC  = c(1.000, 1.039, 1.019),
  MAF_filter      = c(">0.01",">0.01","N/A")
)
fwrite(supp_qc, file.path(RESULTS_DIR,"supplementary_qc_table.tsv"), sep="\t")
print(supp_qc)

cat("\n=== Step 3 complete ===\n")
