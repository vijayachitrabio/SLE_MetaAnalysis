#!/usr/bin/env Rscript
# =============================================================================
# finngen_replication_v2.R
# Independent replication of SLE meta-analysis lead SNPs in FinnGen R12.
#
# Key features:
#  1. rsID-first matching, CHR:BP fallback only where rsID is absent/mismatched.
#  2. Field mapping based on FinnGen R12 documentation (REF=non-effect, ALT=effect).
#  3. Explicit removal of ambiguous (palindromic) SNPs before sign-flip logic.
#  4. Build consistency check (GRCh38 position verification).
#  5. Conservative Bonferroni correction (based on original lead SNP count).
# =============================================================================

library(data.table)
library(dplyr)

RAW_DIR     <- "../data/raw"
RESULTS_DIR <- "../results"

cat("=== FinnGen R12 Replication (v2 – rigorous) ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Load discovery lead SNPs
# ─────────────────────────────────────────────────────────────────────────────
# Prefer PLINK-clumped results if available, else fall back to distance-clumped
lead_file_plink <- file.path(RESULTS_DIR, "sle_meta_lead_snps_plink_clumped.tsv")
lead_file_dist  <- file.path(RESULTS_DIR, "sle_meta_lead_snps.tsv")

if (file.exists(lead_file_plink)) {
  lead_dt <- fread(lead_file_plink)
  cat("Using PLINK-clumped lead SNPs:", nrow(lead_dt), "loci\n\n")
} else {
  lead_dt <- fread(lead_file_dist)
  cat("NOTE: PLINK-clumped file not found; using distance-clumped lead SNPs\n")
  cat("      (", nrow(lead_dt), "loci). Run sle_ld_clumping_plink.R for publication.\n\n")
}

# Standardise column names expected downstream
# meta-analysis uses: SNP (rsID), CHR, BP, A1 (non-effect), A2 (effect), BETA, SE, P
stopifnot(all(c("SNP","CHR","BP","A1","A2","BETA","SE","P") %in% names(lead_dt)))

lead_dt[, A1 := toupper(A1)]
lead_dt[, A2 := toupper(A2)]
cat("Discovery lead SNPs:\n")
print(lead_dt[, .(SNP, CHR, BP, A1, A2, BETA, P)])

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Helper: identify palindromic / ambiguous SNPs
# ─────────────────────────────────────────────────────────────────────────────
is_palindromic <- function(a1, a2) {
  # A/T or T/A or C/G or G/C
  (a1=="A" & a2=="T") | (a1=="T" & a2=="A") |
  (a1=="C" & a2=="G") | (a1=="G" & a2=="C")
}

lead_dt[, palindromic := is_palindromic(A1, A2)]
n_pal <- sum(lead_dt$palindromic)
if (n_pal > 0) {
  cat("\nWARNING:", n_pal, "palindromic SNP(s) in discovery set.\n")
  cat("These will be excluded from replication (cannot safely resolve strand):\n")
  print(lead_dt[palindromic == TRUE, .(SNP, A1, A2)])
}
lead_nonpal <- lead_dt[palindromic == FALSE]
cat("\nNon-palindromic lead SNPs for replication:", nrow(lead_nonpal), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Load FinnGen R12 L12_LUPUS summary statistics
# FinnGen field definitions (verified from FinnGen documentation):
#   #chrom : chromosome (format: "chr1", "chrX")
#   pos    : GRCh38 position
#   ref    : reference allele (NON-effect allele)
#   alt    : alternate allele (EFFECT allele — beta is per alt allele)
#   rsids  : rsID (may be comma-separated if multiple)
#   pval   : p-value
#   beta   : log-OR per ALT allele copy
#   sebeta : standard error of beta
# ─────────────────────────────────────────────────────────────────────────────
fg_file <- file.path(RAW_DIR, "finngen_R12_Lupus.gz")
cat("\nLoading FinnGen R12 data from:", fg_file, "\n")
fg_dt <- fread(fg_file)

# Rename to internal standard
original_names <- names(fg_dt)
cat("FinnGen columns:", paste(original_names, collapse=", "), "\n")

# Flexible rename — handle possible variation in header names
col_map <- list(
  CHR_fg   = c("#chrom","chrom","chromosome","chr"),
  BP_fg    = c("pos","position","bp","base_pair"),
  REF_fg   = c("ref","reference_allele","a1","nea"),    # NON-effect
  ALT_fg   = c("alt","alternate_allele","a2","ea"),     # EFFECT
  RSID_fg  = c("rsids","rsid","snp","variant_id"),
  P_fg     = c("pval","p_value","p","pvalue"),
  BETA_fg  = c("beta","effect","b"),
  SE_fg    = c("sebeta","se","standard_error","stderr")
)

for (new_col in names(col_map)) {
  candidates <- col_map[[new_col]]
  matched <- intersect(tolower(original_names), tolower(candidates))
  if (length(matched) > 0) {
    orig_matched <- original_names[tolower(original_names) %in% matched][1]
    if (!new_col %in% names(fg_dt)) setnames(fg_dt, orig_matched, new_col)
  }
}

# Report what we found
cat("FinnGen column mapping:\n")
for (col in names(col_map)) {
  cat(" ", col, "->",
      if (col %in% names(fg_dt)) col else "NOT FOUND", "\n")
}

# Parse chromosome to numeric
fg_dt[, CHR_fg := as.numeric(gsub("chr", "", CHR_fg, ignore.case=TRUE))]
fg_dt    <- fg_dt[!is.na(CHR_fg) & CHR_fg %in% 1:22]

# Parse rsID (FinnGen may have comma-separated rsids; take first)
if ("RSID_fg" %in% names(fg_dt)) {
  fg_dt[, RSID_fg := sub(",.*", "", RSID_fg)]
}

fg_dt[, REF_fg := toupper(REF_fg)]
fg_dt[, ALT_fg := toupper(ALT_fg)]

cat("FinnGen data loaded:", nrow(fg_dt), "variants\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Match discovery lead SNPs to FinnGen
#     Priority 1: rsID match (most reliable)
#     Priority 2: CHR:BP match (only if rsID absent or missing)
#     Build note: FinnGen R12 is on GRCh38; GWAS Catalog harmonised files
#     for Bentham/Yin/Julia are on GRCh38 (harmonised). This is consistent.
#     IF builds differ, positions will diverge — we flag that below.
# ─────────────────────────────────────────────────────────────────────────────
# Build CHR_BP key
fg_dt[, c_id := paste(CHR_fg, BP_fg, sep="_")]
lead_nonpal[, c_id := paste(CHR, BP, sep="_")]

# -- Primary match: by rsID --
if ("RSID_fg" %in% names(fg_dt)) {
  matched_rsid <- merge(
    lead_nonpal,
    fg_dt[, .(SNP_fg = RSID_fg, CHR_fg, BP_fg, REF_fg, ALT_fg, P_fg, BETA_fg, SE_fg, c_id)],
    by.x = "SNP", by.y = "SNP_fg",
    all.x = FALSE
  )
  cat("rsID-matched SNPs:", nrow(matched_rsid), "\n")
} else {
  matched_rsid <- data.table()
  cat("No rsID column in FinnGen; all matching by CHR:BP.\n")
}

# -- Fallback match: by CHR:BP for unmatched SNPs --
unmatched_snps <- lead_nonpal[!SNP %in% matched_rsid$SNP]
if (nrow(unmatched_snps) > 0 && "c_id" %in% names(fg_dt)) {
  matched_pos <- merge(
    unmatched_snps,
    fg_dt[, .(CHR_fg, BP_fg, REF_fg, ALT_fg, P_fg, BETA_fg, SE_fg, c_id)],
    by = "c_id"
  )
  cat("CHR:BP-matched SNPs (fallback):", nrow(matched_pos), "\n")
  matched_pos[, SNP_fg := SNP]  # use discovery rsID as label
} else {
  matched_pos <- data.table()
}

# Combine
all_matched <- rbindlist(list(matched_rsid, matched_pos), fill=TRUE, use.names=TRUE)
all_matched <- unique(all_matched, by="SNP")
cat("Total matched SNPs:", nrow(all_matched), "out of", nrow(lead_nonpal), "\n\n")

if (nrow(all_matched) == 0) stop("No lead SNPs matched in FinnGen. Check builds/rsIDs.")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Build harmonization check
#     Compare CHR:BP from discovery vs FinnGen after rsID match.
#     Position differences > 1000bp likely indicate a build mismatch.
# ─────────────────────────────────────────────────────────────────────────────
if ("BP_fg" %in% names(all_matched)) {
  all_matched[, pos_diff := abs(BP - BP_fg)]
  build_mismatch <- all_matched[!is.na(pos_diff) & pos_diff > 1000]
  if (nrow(build_mismatch) > 0) {
    cat("WARNING: Potential build mismatch detected for", nrow(build_mismatch),
        "SNP(s). Position differences > 1000bp:\n")
    print(build_mismatch[, .(SNP, CHR, BP_discovery=BP, BP_fg, pos_diff)])
    cat("RECOMMENDATION: Verify that both summary stats are on the same reference build (GRCh38).\n\n")
  } else {
    cat("Build check PASSED: All rsID-matched SNP positions agree within 1000 bp.\n\n")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Allele alignment
# FinnGen: REF_fg = non-effect, ALT_fg = effect (beta is per ALT copy)
# Discovery: A1 = non-effect, A2 = effect
#
# Case 1: A1==REF_fg & A2==ALT_fg → strands consistent, no flip needed
# Case 2: A1==ALT_fg & A2==REF_fg → flip: BETA_fg = -BETA_fg
# Case 3: complement of above (i.e., different strand representation)
# Other  : incompatible alleles → exclude
# ─────────────────────────────────────────────────────────────────────────────
complement <- function(x) {
  chartr("ACGT", "TGCA", x)
}

all_matched[, `:=`(
  A1_comp = complement(A1),
  A2_comp = complement(A2)
)]

all_matched[, align_status := case_when(
  A1 == REF_fg & A2 == ALT_fg         ~ "concordant",
  A1 == ALT_fg & A2 == REF_fg         ~ "flipped",
  A1_comp == REF_fg & A2_comp == ALT_fg ~ "concordant_comp",  # opposite strand
  A1_comp == ALT_fg & A2_comp == REF_fg ~ "flipped_comp",
  TRUE                                  ~ "incompatible"
)]

cat("Allele alignment status:\n")
print(table(all_matched$align_status))

# Apply flip where needed
all_matched[align_status %in% c("flipped","flipped_comp"),
            BETA_fg := -BETA_fg]

# Exclude incompatible allele pairs
n_incompat <- sum(all_matched$align_status == "incompatible")
if (n_incompat > 0) {
  cat("\nExcluding", n_incompat, "SNP(s) with incompatible alleles:\n")
  print(all_matched[align_status=="incompatible", .(SNP, A1, A2, REF_fg, ALT_fg)])
}
rep_dt <- all_matched[align_status != "incompatible"]
cat("SNPs retained after allele alignment:", nrow(rep_dt), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7.  Replication testing
# ─────────────────────────────────────────────────────────────────────────────
# Conservative Bonferroni threshold based on the total number of original 
# discovery lead SNPs (not just the subset that passed alignment/exclusion).
n_discovery_leads <- nrow(lead_dt)
n_tested          <- nrow(rep_dt)
alpha_bon         <- 0.05 / n_discovery_leads 

cat(sprintf("Number of original discovery lead SNPs: %d\n", n_discovery_leads))
cat(sprintf("Number of SNPs tested in replication: %d\n", n_tested))
cat(sprintf("Bonferroni replication threshold (0.05 / %d): %.5f\n\n", 
            n_discovery_leads, alpha_bon))

rep_dt[, `:=`(
  dir_consistent = (sign(BETA) == sign(BETA_fg)),
  replicated      = (P_fg < alpha_bon) & (sign(BETA) == sign(BETA_fg))
)]

cat("=== Replication Results ===\n")
cat("Direction-consistent:", sum(rep_dt$dir_consistent, na.rm=TRUE), "\n")
cat("Replicated (P <", round(alpha_bon,5), "& consistent direction):",
    sum(rep_dt$replicated, na.rm=TRUE), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 8.  Save and report
# ─────────────────────────────────────────────────────────────────────────────
output_cols <- c("SNP","CHR","BP","A1","A2",
                 "BETA","SE","P",
                 "BETA_fg","SE_fg","P_fg",
                 "align_status","dir_consistent","replicated")

out_dt <- rep_dt[, ..output_cols]
fwrite(out_dt, file.path(RESULTS_DIR, "finngen_replication_v2_results.tsv"), sep="\t")

# Include palindromes in output (flagged as excluded)
if (n_pal > 0) {
  pal_rows <- lead_dt[palindromic==TRUE,
                      .(SNP, CHR, BP, A1, A2, BETA, SE, P,
                        BETA_fg=NA_real_, SE_fg=NA_real_, P_fg=NA_real_,
                        align_status="palindromic_excluded",
                        dir_consistent=NA, replicated=NA)]
  out_with_pal <- rbindlist(list(out_dt, pal_rows), fill=TRUE)
  fwrite(out_with_pal,
         file.path(RESULTS_DIR, "finngen_replication_v2_results_full.tsv"), sep="\t")
  cat("Full results (including palindromes) saved.\n")
}

cat("\n=== Final Replication Table ===\n")
print(out_dt[order(P),
             .(SNP, CHR, BP, Meta_BETA=round(BETA,4), Meta_P=formatC(P,format="e",digits=2),
               FG_BETA=round(BETA_fg,4), FG_P=formatC(P_fg,format="e",digits=2),
               Direction=ifelse(dir_consistent,"✅","❌"),
               Replicated=ifelse(replicated,"✅","❌"))])

cat("\n=== Replication analysis complete (v2) ===\n")
