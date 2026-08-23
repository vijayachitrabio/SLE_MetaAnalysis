#!/usr/bin/env Rscript
# revisions_todays_data/run_lava_bentham_only_sensitivity.R
#
# ZERO-SAMPLE-OVERLAP SENSITIVITY ANALYSIS for the bivariate LAVA cross-trait
# results (Reviewer 1, Point 5).
#
# Rationale: the primary bivariate analysis used the full SLE discovery
# meta-analysis (Bentham 2015 + FinnGen R12). All three comparator traits
# (RA, SSc, Sjogren) are FinnGen R12-derived and therefore share control
# samples with the FinnGen component of the SLE discovery set. Unmodelled
# sample overlap biases LAVA local genetic correlations upward.
#
# This script repeats the bivariate analysis using Bentham 2015 ALONE as the
# SLE input. Bentham is a UK/European non-Finnish cohort with ZERO sample
# overlap with FinnGen R12 by construction, so no overlap modelling is needed.
#
# IMPORTANT CAVEAT (reported alongside results): Bentham alone has far fewer
# controls (15,991 vs 354,277), so effective N drops from ~32,900 to ~19,900
# (~40%). Fewer loci will pass LAVA's univariate h2 threshold. This test can
# CONFIRM robustness convincingly, but an attenuated result is ambiguous
# between "overlap was inflating it" and "we lost power".
#
# ALLELE ALIGNMENT NOTE: the Bentham GWAS Catalog file carries both raw
# (effect_allele/beta) and harmonized (hm_effect_allele/hm_beta) column sets.
# These differ for variants that were strand/allele-flipped during
# harmonization (e.g. rs4607995: hm_effect=T/beta=+0.030 vs effect=C/beta=-0.030).
# We therefore use the hm_* set CONSISTENTLY (hm_rsid, hm_effect_allele,
# hm_other_allele, hm_beta) and never mix it with the raw columns.
# standard_error is allele-flip invariant so the raw column is safe.

suppressPackageStartupMessages({
    library(data.table)
    library(LAVA)
})

setDTthreads(0)

OUT_DIR <- "revisions_todays_data"
MUNGED  <- "data/munged"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BENTHAM_RAW <- "data/raw/Bentham_2015_SLE.h.tsv.gz"
BENTHAM_OUT <- file.path(MUNGED, "SLE_Bentham_only_lava_input.txt")

# Bentham et al. 2015 European SLE GWAS
N_CASES    <- 7219
N_CONTROLS <- 15991
N_TOTAL    <- N_CASES + N_CONTROLS   # 23,210

cat("=========================================\n")
cat("Bentham-only zero-overlap LAVA sensitivity\n")
cat("=========================================\n\n")

## ---------------------------------------------------------------- 1. Munge
if (!file.exists(BENTHAM_OUT)) {
    cat("Munging Bentham-only SLE sumstats (this is the slow step)...\n")
    b <- fread(
        BENTHAM_RAW,
        select = c("hm_rsid", "hm_effect_allele", "hm_other_allele",
                   "hm_beta", "standard_error", "p_value",
                   "effect_allele_frequency"),
        showProgress = TRUE
    )
    cat(sprintf("  Read %s rows.\n", format(nrow(b), big.mark = ",")))

    # Keep only rows where the harmonized set is complete
    b <- b[!is.na(hm_rsid) & hm_rsid != "" &
           !is.na(hm_beta) & !is.na(standard_error) & !is.na(p_value) &
           standard_error > 0]
    cat(sprintf("  %s rows with complete harmonized stats.\n",
                format(nrow(b), big.mark = ",")))

    # MAF filter (AF is allele-flip invariant for a symmetric MAF threshold)
    b <- b[is.na(effect_allele_frequency) |
           (effect_allele_frequency > 0.01 & effect_allele_frequency < 0.99)]
    cat(sprintf("  %s rows after MAF filter.\n", format(nrow(b), big.mark = ",")))

    # LAVA format: A1 = effect allele, Z = beta / se
    b[, Z := hm_beta / standard_error]
    b[, N := N_TOTAL]
    setnames(b,
             c("hm_rsid", "hm_effect_allele", "hm_other_allele", "p_value"),
             c("SNP", "A1", "A2", "P"))

    # Same deduplication policy as harmonize_lava_traits.R:
    # drop allele-conflicting duplicate groups, else keep lowest P.
    dup_snps <- unique(b$SNP[duplicated(b$SNP)])
    if (length(dup_snps) > 0) {
        cat(sprintf("  Resolving %s duplicated RSIDs...\n",
                    format(length(dup_snps), big.mark = ",")))
        b_uniq <- b[!SNP %in% dup_snps]
        b_dup  <- b[SNP %in% dup_snps]
        b_dup_res <- b_dup[, {
            if (uniqueN(paste(A1, A2, sep = "_")) > 1) .SD[0] else .SD[which.min(P)]
        }, by = SNP]
        b <- rbind(b_uniq[, .(SNP, A1, A2, Z, N, P)],
                   b_dup_res[, .(SNP, A1, A2, Z, N, P)])
    } else {
        b <- b[, .(SNP, A1, A2, Z, N, P)]
    }

    stopifnot(!any(duplicated(b$SNP)))
    fwrite(b, BENTHAM_OUT, sep = "\t")
    cat(sprintf("  Saved %s unique variants -> %s\n\n",
                format(nrow(b), big.mark = ","), BENTHAM_OUT))
    rm(b); gc()
} else {
    cat("Reusing existing munged Bentham-only input.\n\n")
}

## ------------------------------------------------------ 2. LAVA input file
info_file <- file.path(MUNGED, "lava_bentham_sensitivity_info.txt")
info <- data.table(
    phenotype = c("SLE_Bentham", "RA", "SSc", "Sjogren"),
    cases     = c(N_CASES, 16314, 854, 3309),
    controls  = c(N_CONTROLS, 298801, 483406, 480951),
    filename  = c(BENTHAM_OUT,
                  file.path(MUNGED, "RA_lava_input.txt"),
                  file.path(MUNGED, "SSc_lava_input.txt"),
                  file.path(MUNGED, "Sjogren_lava_input.txt"))
)
fwrite(info, info_file, sep = "\t")
cat("Wrote LAVA info file:\n"); print(info); cat("\n")

## --------------------------------------------------------- 3. LAVA process
REF_PREFIX <- "reference_data/g1000_eur"
BLOCK_FILE <- "reference_data/LAVA_s2500_m25_f1_w200.blocks"
LOCI_FILE  <- "results_extracted/final_lava_consolidated_loci.tsv"

cat("Processing LAVA input (SNP alignment)...\n")
input <- process.input(
    input.info.file     = info_file,
    sample.overlap.file = NULL,   # legitimately NULL: zero overlap by construction
    ref.prefix          = REF_PREFIX,
    phenos              = c("SLE_Bentham", "RA", "SSc", "Sjogren")
)

# NOTE: the .blocks file has lowercase headers (chr/start/stop) and NO LOC
# column. The primary bivariate and conditional scripts both derive LOC as the
# sequential ROW INDEX of the blocks file. We replicate that convention exactly
# so LOC ids correspond 1:1 with the primary run; using LAVA::read.loci() here
# would fail and any other indexing would silently misalign loci.
loci_all <- fread(BLOCK_FILE)
setnames(loci_all, names(loci_all), toupper(names(loci_all)))
if (!"LOC" %in% names(loci_all)) loci_all$LOC <- seq_len(nrow(loci_all))

target   <- fread(LOCI_FILE)
target_locs <- as.numeric(target$LAVA_LOC)
target_locs <- sort(unique(target_locs[!is.na(target_locs) & target_locs != ""]))
target_locs <- setdiff(target_locs, 961)  # exclude chr6p21.3/MHC CLIC1 block
cat(sprintf("\nTargeting %d prioritized non-MHC loci (LOC 961/MHC excluded).\n\n",
            length(target_locs)))

## ------------------------------------------------------------ 4. Bivariate
bivar_res <- list()
univ_res  <- list()
skipped   <- list()

for (i in seq_along(target_locs)) {
    loc_id <- target_locs[i]
    cat(sprintf("   [%d/%d] LOC %d...\n", i, length(target_locs), loc_id))

    row <- loci_all[LOC == loc_id]
    if (nrow(row) == 0) {
        skipped[[length(skipped) + 1]] <- data.table(LOC = loc_id,
            Reason = "Locus not found in LAVA block file")
        next
    }

    locus <- tryCatch(process.locus(row, input), error = function(e) NULL)
    if (is.null(locus)) {
        skipped[[length(skipped) + 1]] <- data.table(LOC = loc_id,
            Reason = "process.locus failed or no PCs retained")
        next
    }

    # Match the primary run: no `target` arg, then filter to SLE pairs below.
    out <- tryCatch(run.univ.bivar(locus), error = function(e) NULL)
    if (is.null(out)) {
        skipped[[length(skipped) + 1]] <- data.table(LOC = loc_id,
            Reason = "run.univ.bivar failed")
        next
    }

    if (!is.null(out$univ) && nrow(out$univ) > 0) {
        u <- as.data.table(out$univ); u[, LOC := loc_id]
        univ_res[[length(univ_res) + 1]] <- u
    }

    # LAVA returns the bivariate table as $biv (the primary script uses $biv).
    bvt <- if (!is.null(out$biv)) out$biv else out$bivar
    if (!is.null(bvt) && nrow(bvt) > 0) {
        b <- as.data.table(bvt); b[, LOC := loc_id]
        bivar_res[[length(bivar_res) + 1]] <- b
    } else {
        avail <- if (!is.null(out$univ)) paste(out$univ$phen, collapse = ",") else ""
        skipped[[length(skipped) + 1]] <- data.table(LOC = loc_id,
            Reason = paste0("No bivariate result (univ-eligible: ", avail, ")"))
    }
}

## --------------------------------------------------------------- 5. Output
if (length(univ_res) > 0) {
    fwrite(rbindlist(univ_res, fill = TRUE),
           file.path(OUT_DIR, "lava_bentham_sensitivity_univ.csv"))
    cat(sprintf("\nSaved %d univariate rows.\n",
                nrow(rbindlist(univ_res, fill = TRUE))))
}
if (length(skipped) > 0) {
    fwrite(rbindlist(skipped, fill = TRUE),
           file.path(OUT_DIR, "lava_bentham_sensitivity_skip_log.tsv"), sep = "\t")
    cat(sprintf("Logged %d skipped loci.\n", nrow(rbindlist(skipped, fill = TRUE))))
}

if (length(bivar_res) == 0) {
    cat("\n[WARN] No bivariate results produced. Check the skip log.\n")
    quit(status = 0)
}

bv <- rbindlist(bivar_res, fill = TRUE)
bv <- bv[phen1 == "SLE_Bentham" | phen2 == "SLE_Bentham"]
bv[, Secondary_Trait := ifelse(phen1 == "SLE_Bentham", phen2, phen1)]
bv[, FDR := p.adjust(p, method = "BH")]

ann <- unique(target[, .(LOC = as.integer(LAVA_LOC), Gene, Region, RSID)])
bv <- merge(bv, ann, by = "LOC", all.x = TRUE)

setorder(bv, p)
fwrite(bv, file.path(OUT_DIR, "lava_bentham_sensitivity_bivariate.csv"))

cat(sprintf("\n[SUCCESS] %d bivariate tests (Bentham-only SLE).\n", nrow(bv)))
cat(sprintf("          %d significant at FDR < 0.05.\n", sum(bv$FDR < 0.05)))
cat(sprintf("          Saved -> %s\n",
            file.path(OUT_DIR, "lava_bentham_sensitivity_bivariate.csv")))

## ------------------------------------------- 6. Side-by-side vs primary run
primary_f <- "results/lava_crosstrait_results_sle_only.tsv"
if (file.exists(primary_f)) {
    prim <- fread(primary_f)
    prim[, LOC := as.integer(LOC)]
    prim_s <- prim[, .(LOC, Secondary_Trait,
                       rho_primary = Genetic_Correlation_Rho,
                       FDR_primary = FDR)]
    sens_s <- bv[, .(LOC, Secondary_Trait,
                     rho_bentham = round(rho, 4),
                     FDR_bentham = FDR, Gene, RSID)]
    cmp <- merge(prim_s, sens_s, by = c("LOC", "Secondary_Trait"), all = TRUE)
    cmp[, rho_delta := round(rho_bentham - rho_primary, 4)]
    cmp[, testable_in_both := !is.na(rho_primary) & !is.na(rho_bentham)]
    cmp[, sig_both := !is.na(FDR_primary) & !is.na(FDR_bentham) &
                      FDR_primary < 0.05 & FDR_bentham < 0.05]
    setorder(cmp, -testable_in_both, FDR_bentham, na.last = TRUE)
    fwrite(cmp, file.path(OUT_DIR, "lava_overlap_sensitivity_comparison.tsv"),
           sep = "\t")

    nb <- sum(cmp$testable_in_both)
    cat("\n--- Primary vs Bentham-only comparison ---\n")
    cat(sprintf("  Tests in primary run           : %d\n", nrow(prim_s)))
    cat(sprintf("  Tests in Bentham-only run      : %d\n", nrow(sens_s)))
    cat(sprintf("  Testable in BOTH               : %d\n", nb))
    if (nb > 0) {
        both <- cmp[testable_in_both == TRUE]
        cat(sprintf("  FDR-sig in both                : %d\n", sum(both$sig_both)))
        cat(sprintf("  Mean rho (primary)             : %.3f\n",
                    mean(both$rho_primary, na.rm = TRUE)))
        cat(sprintf("  Mean rho (Bentham-only)        : %.3f\n",
                    mean(both$rho_bentham, na.rm = TRUE)))
        cat(sprintf("  Mean change in rho             : %+.3f\n",
                    mean(both$rho_delta, na.rm = TRUE)))
        cat(sprintf("  Directionally consistent       : %d / %d\n",
                    sum(sign(both$rho_primary) == sign(both$rho_bentham),
                        na.rm = TRUE), nb))
    }
    cat(sprintf("  Saved -> %s\n",
                file.path(OUT_DIR, "lava_overlap_sensitivity_comparison.tsv")))
}

cat("\nDone.\n")
