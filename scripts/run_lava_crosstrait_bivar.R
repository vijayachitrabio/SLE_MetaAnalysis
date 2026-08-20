#!/usr/bin/env Rscript
# scripts/run_lava_crosstrait_bivar.R
# LAVA Bivariate Analysis for SLE vs RA, SSc, Sjögren

suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
})

cat("=========================================\n")
cat("Starting LAVA Cross-Trait Analysis\n")
cat("=========================================\n")

# --- Configuration ---
WD <- "."
setwd(WD)
REF_PREFIX   <- "reference_data/g1000_eur" # Will need to ensure this points to valid bed/bim/fam
BLOCKS_FILE  <- "reference_data/LAVA_s2500_m25_f1_w200.blocks"

# Prioritized Loci list (47 genome-wide significant SLE loci)
LAVA_LOCI_FILE <- "results_extracted/final_lava_consolidated_loci.tsv"
if (file.exists(LAVA_LOCI_FILE)) {
    sle_loci <- fread(LAVA_LOCI_FILE)
    target_loc_ids <- unique(sle_loci$LAVA_LOC)
    # Filter out empty or missing LOC IDs
    target_loc_ids <- as.numeric(target_loc_ids[!is.na(target_loc_ids) & target_loc_ids != ""])
} else {
    cat("[WARNING] Could not find final_lava_consolidated_loci.tsv. Will run all blocks.\n")
    target_loc_ids <- NULL
}

# ============================================================
# 0. Create SLE Discovery Input
# ============================================================
cat("\n--- Formatting Trait 1: SLE Discovery Meta ---\n")
d1 <- fread("results/discovery_meta_results.tsv",
            select = c("CHR", "BP", "EA", "OA", "RSID", "P_meta", "BETA_meta", "SE_meta"))
setnames(d1, c("CHR", "BP", "A1", "A2", "SNP", "P", "BETA", "SE"))
d1 <- d1[!is.na(P) & !is.na(BETA) & !is.na(SE)]
d1 <- d1[P > 0 & P <= 1 & SE > 0 & is.finite(BETA) & is.finite(SE)]
d1[, Z := BETA / SE]
d1[, N := 362694] # N_DISCO (8417 cases + 354277 controls)
fwrite(d1[, .(SNP, A1, A2, Z, P, N)], "data/munged/SLE_lava_input.txt", sep="\t")
rm(d1); gc()

# ============================================================
# 1. Create LAVA Input Info File
# ============================================================
# Note: Sample sizes for SLE are from discovery meta (388,655).
# RA (315,115), SSc (484,260), Sjogren (484,260)
info_content <- data.frame(
  phenotype = c("SLE_Discovery", "RA", "SSc", "Sjogren"),
  cases     = c(8417, 16314, 854, 3309), 
  controls  = c(354277, 298801, 483406, 480951),  
  filename  = c("data/munged/SLE_lava_input.txt",
                "data/munged/RA_lava_input.txt",
                "data/munged/SSc_lava_input.txt",
                "data/munged/Sjogren_lava_input.txt")
)
cat("\n[INFO] LAVA cross-trait input info:\n")
print(info_content)
write.table(info_content, "data/munged/lava_crosstrait_info.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================
# 2. Process Input with LAVA
# ============================================================
cat("\nProcessing LAVA cross-trait input (SNP alignment)...\n")
input <- process.input(
  input.info.file    = "data/munged/lava_crosstrait_info.txt",
  sample.overlap.file = NULL,
  ref.prefix         = REF_PREFIX
)

# ============================================================
# 3. Load Blocks and Run Analysis Loop
# ============================================================
loci <- fread(BLOCKS_FILE)
setnames(loci, names(loci), toupper(names(loci)))
if (!"LOC" %in% names(loci)) loci$LOC <- seq_len(nrow(loci))

if (!is.null(target_loc_ids)) {
    loci <- loci[LOC %in% target_loc_ids]
}

total_loci <- nrow(loci)
cat(sprintf("Starting modeling loop for %d prioritized loci...\n", total_loci))

results_list <- list()
univ_results_list <- list()

for (i in seq_len(total_loci)) {
  loc_id <- loci$LOC[i]
  cat(sprintf("   [%d/%d] Processing LOC %s...\n", i, total_loci, loc_id))
  
  locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
  if (is.null(locus)) next
  
  res <- tryCatch(run.univ.bivar(locus), error = function(e) NULL)
  
  if (!is.null(res) && !is.null(res$univ)) {
    univ_df <- as.data.table(res$univ)
    univ_df[, LOC := loc_id]
    univ_results_list[[i]] <- univ_df
  }
  
  if (!is.null(res) && !is.null(res$biv)) {
    res_df <- as.data.table(res$biv)
    res_df[, LOC := loc_id]
    results_list[[i]] <- res_df
  }
}

if (length(univ_results_list) > 0) {
  final_univ <- rbindlist(univ_results_list)
  fwrite(final_univ, "results/lava_univ_results.csv")
}


# ============================================================
# 4. Final Results Consolidation
# ============================================================
if (length(results_list) > 0) {
  final_results <- rbindlist(results_list)
  fwrite(final_results, "results/lava_crosstrait_bivariate.csv")
  cat(sprintf("\n[SUCCESS] %d bivariate results saved to: results/lava_crosstrait_bivariate.csv\n",
              nrow(final_results)))
} else {
  cat("[ERROR] No bivariate results were produced. Check LD reference and input data.\n")
}

cat("Done!\n")
