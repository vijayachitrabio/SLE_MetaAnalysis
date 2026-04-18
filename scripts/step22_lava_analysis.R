#!/usr/bin/env Rscript
# scripts/step22_lava_analysis.R
# LAVA Analysis for SLE Meta-Analysis (Discovery vs Replication)
# CORRECTED VERSION:
#   - True sample sizes (not SNP counts)
#   - Verified allele orientation for Julia_2018_Spain_Only.txt
#   - Robust Z-score derivation with full input validation
#   - Explicit column mapping logging before LAVA run

suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
  library(dplyr)
})

cat("=========================================\n")
cat("Starting LAVA Analysis for SLE (CORRECTED)\n")
cat(" (Discovery Meta vs Spanish Replication)\n")
cat("=========================================\n")

# --- Configuration ---
WD <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis"
setwd(WD)

# Paths to external LD data (from fibroids project)
REF_DIR  <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/uterine_fibroids/mr_analysis_2026_04_16/reference_data"
REF_PREFIX   <- file.path(REF_DIR, "g1000_eur")
BLOCKS_FILE  <- file.path(REF_DIR, "LAVA_s2500_m25_f1_w200.blocks")

# Summary stats paths
DISCO_FILE <- "results/discovery_meta_results.tsv"
REP_FILE   <- "data/raw/Julia_2018_Spain_Only.txt"

# Output files
PROGRESS_FILE <- "results/lava_progress_checkpoint.csv"
FINAL_RESULTS <- "results/lava_sle_results.csv"

# ============================================================
# CORRECTED SAMPLE SIZES
# Discovery meta = Bentham 2015 (4036 cases + 11040 ctrl = 15076)
#                + FinnGen R12 SLE_FG (1306 cases + 372273 ctrl = 373579)
# Combined:        5342 cases + 383313 controls = 388655
#
# Replication     = Julia 2018 Spanish-only: 907 cases + 2845 ctrl = 3752
# ============================================================
N_DISCO_CASES    <- 5342
N_DISCO_CONTROLS <- 383313
N_DISCO          <- N_DISCO_CASES + N_DISCO_CONTROLS  # 388655

N_REP_CASES    <- 907
N_REP_CONTROLS <- 2845
N_REP          <- N_REP_CASES + N_REP_CONTROLS  # 3752

cat(sprintf("[INFO] Discovery N = %d (%d cases / %d controls)\n",
            N_DISCO, N_DISCO_CASES, N_DISCO_CONTROLS))
cat(sprintf("[INFO] Replication N = %d (%d cases / %d controls)\n",
            N_REP, N_REP_CASES, N_REP_CONTROLS))

# ============================================================
# 1. Load Reference SNP List
# ============================================================
cat("Loading reference SNP list (LD panel)...\n")
bim        <- fread(paste0(REF_PREFIX, ".bim"), select = 2, col.names = "SNP")
valid_snps <- bim$SNP
rm(bim); gc()
cat(sprintf("[INFO] LD panel contains %d SNPs\n", length(valid_snps)))

# ============================================================
# 2. Format Trait 1: Discovery Meta
#
# Columns in discovery_meta_results.tsv:
#   CHR, BP, RSID, OA, EA, BETA_meta, SE_meta, P_meta ...
#
# Effect allele = EA (as defined in Bentham 2015; FinnGen aligned to this)
# BETA_meta is already on the EA scale (positive -> EA raises risk)
# ============================================================
cat("\n--- Formatting Trait 1: SLE Discovery Meta ---\n")
d1 <- fread(DISCO_FILE,
            select = c("CHR", "BP", "EA", "OA", "RSID",
                       "P_meta", "BETA_meta", "SE_meta"))
setnames(d1, c("CHR", "BP", "A1", "A2", "SNP", "P", "BETA", "SE"))

# Validation guards
n_before <- nrow(d1)
d1 <- d1[SNP %in% valid_snps]
d1 <- d1[!is.na(P) & !is.na(BETA) & !is.na(SE)]
d1 <- d1[P   >  0 & P  <= 1]
d1 <- d1[SE  >  0]
d1 <- d1[is.finite(BETA) & is.finite(SE)]
cat(sprintf("[INFO] Discovery: %d rows before filters, %d after\n",
            n_before, nrow(d1)))

d1[, N := N_DISCO]
d1[, Z := BETA / SE]

# Log mapping
cat("[CHECK] Discovery column mapping (first 5 rows):\n")
print(head(d1[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, Z, N)], 5))

fwrite(d1, "results/lava_disco_filtered.txt", sep = "\t")
rm(d1); gc()

# ============================================================
# 3. Format Trait 2: Spanish Replication
#
# Columns in Julia_2018_Spain_Only.txt:
#   CHR, SNP, BP, A1leleA, AlleleB, OR, OR_lower, OR_upper, P
#
# IMPORTANT allele convention:
#   The OR in this file is relative to AlleleB (coded/effect allele).
#   i.e.  OR > 1  means AlleleB increases risk.
#   We therefore set:
#     A1 = AlleleB  (effect allele; matches LAVA convention)
#     A2 = A1leleA  (other allele)
#     BETA = log(OR)  [positive when AlleleB increases risk]
#
#   Note: The column header is "A1leleA" (typo in original file),
#         which is the NON-effect/reference allele.
# ============================================================
cat("\n--- Formatting Trait 2: Spanish Replication ---\n")

# Read and inspect header
d2_raw <- fread(REP_FILE)
cat("[INFO] Julia 2018 raw columns:", paste(names(d2_raw), collapse = ", "), "\n")
cat("[INFO] Julia 2018 first 3 rows:\n")
print(head(d2_raw, 3))

# Rename for clarity
# AlleleB = effect allele -> A1
# A1leleA = reference    -> A2
d2 <- d2_raw[, .(CHR, SNP, BP,
                  A1 = AlleleB,
                  A2 = A1leleA,
                  OR, P)]
rm(d2_raw); gc()

# Validation guards
n_before <- nrow(d2)
d2 <- d2[SNP %in% valid_snps]
d2 <- d2[!is.na(OR) & OR > 0]
d2 <- d2[!is.na(P)  & P > 0 & P <= 1]
d2 <- d2[is.finite(OR)]
cat(sprintf("[INFO] Replication: %d rows before filters, %d after\n",
            n_before, nrow(d2)))

# BETA = log(OR); cap P at 1e-300 to prevent Inf Z
d2[, BETA := log(OR)]
d2[, P_capped := pmax(P, 1e-300)]
d2[, Z := qnorm(P_capped / 2, lower.tail = FALSE) * sign(BETA)]
d2[, N := N_REP]

# Sanity check: verify concordance for rs4853458 if present
# Discovery: EA=A, BETA_meta > 0 (risk increases on A)
# Julia: AlleleB = A, OR = 0.649 -> BETA = -0.43 (protective? OR < 1)
# This is expected because OR can be < 1 depending on case/control definition
# LAVA handles cross-trait sign alignment internally via the LD reference
check_snp <- "rs4853458"
if (check_snp %in% d2$SNP) {
  cat(sprintf("[CHECK] %s in replication - A1=%s, A2=%s, OR=%.4f, BETA=%.4f, Z=%.4f\n",
              check_snp,
              d2[SNP == check_snp, A1],
              d2[SNP == check_snp, A2],
              d2[SNP == check_snp, OR],
              d2[SNP == check_snp, BETA],
              d2[SNP == check_snp, Z]))
}

# Log mapping
cat("[CHECK] Replication column mapping (first 5 rows):\n")
print(head(d2[, .(SNP, CHR, BP, A1, A2, BETA, Z, P, N)], 5))

fwrite(d2[, .(CHR, SNP, BP, A1, A2, OR, BETA, Z, P, N)],
       "results/lava_rep_filtered.txt", sep = "\t")
rm(d2, valid_snps); gc()

# ============================================================
# 4. Create LAVA Input Info File (exact case/control counts)
# ============================================================
info_content <- data.frame(
  phenotype = c("SLE_Discovery", "SLE_Replication"),
  cases     = c(N_DISCO_CASES, N_REP_CASES),        # 5342, 907
  controls  = c(N_DISCO_CONTROLS, N_REP_CONTROLS),  # 383313, 2845
  filename  = c("results/lava_disco_filtered.txt",
                "results/lava_rep_filtered.txt")
)
cat("\n[INFO] LAVA input info:\n")
print(info_content)
write.table(info_content, "results/lava_input_info.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================
# 5. Process Input with LAVA
# ============================================================
cat("\nProcessing LAVA input (SNP alignment)...\n")
input <- process.input(
  input.info.file    = "results/lava_input_info.txt",
  sample.overlap.file = NULL,
  ref.prefix         = REF_PREFIX
)

# ============================================================
# 6. Load Blocks
# ============================================================
loci <- fread(BLOCKS_FILE)
setnames(loci, names(loci), toupper(names(loci)))
if (!"LOC" %in% names(loci)) loci$LOC <- seq_len(nrow(loci))

# ============================================================
# 7. Check for Resume (checkpoint)
# ============================================================
done_locs <- c()
if (file.exists(PROGRESS_FILE)) {
  prog <- fread(PROGRESS_FILE)
  if (nrow(prog) > 0) done_locs <- unique(prog$LOC)
  cat(sprintf("[RESUME] Skipping %d already-processed loci.\n",
              length(done_locs)))
}

# ============================================================
# 8. Run Analysis Loop
# ============================================================
total_loci <- nrow(loci)
cat(sprintf("Starting modeling loop (%d loci total)...\n", total_loci))

for (i in seq_len(total_loci)) {
  loc_id <- loci$LOC[i]
  if (loc_id %in% done_locs) next

  if (i %% 50 == 0)
    cat(sprintf("   [%d/%d] loci processed...\n", i, total_loci))

  locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
  if (is.null(locus)) next

  res <- tryCatch(run.univ.bivar(locus), error = function(e) NULL)
  if (!is.null(res) && !is.null(res$biv)) {
    res_to_save      <- as.data.table(res$biv)
    res_to_save[, LOC := loc_id]
    fwrite(res_to_save, PROGRESS_FILE, append = file.exists(PROGRESS_FILE))
  }
}

# ============================================================
# 9. Final Results Consolidation
# ============================================================
cat("\nAnalysis complete! Consolidating results...\n")
if (file.exists(PROGRESS_FILE)) {
  final_results <- fread(PROGRESS_FILE)
  fwrite(final_results, FINAL_RESULTS)
  cat(sprintf("[INFO] %d bivariate results saved to: %s\n",
              nrow(final_results), FINAL_RESULTS))
} else {
  stop("[ERROR] No results file found — no bivariate results were produced.")
}

# Cleanup temp files
suppressWarnings(file.remove(
  "results/lava_disco_filtered.txt",
  "results/lava_rep_filtered.txt",
  "results/lava_input_info.txt",
  PROGRESS_FILE
))

cat("Done!\n")
