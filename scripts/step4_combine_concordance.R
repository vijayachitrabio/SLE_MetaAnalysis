#!/usr/bin/env Rscript
# step4_combine_concordance.R
# Combine: meta-analysis lead SNPs + FinnGen replication
# into a single European-ancestry replication concordance table for manuscript Table 2.
# Cross-ancestry look-ups (Yin 2022, Sakaue 2021) excluded — European-only analysis.

library(data.table)

RESULTS_DIR <- "../results"

cat("=== Step 4: European Replication Concordance Table ===\n\n")

# ── Load all results ──────────────────────────────────────────────────────────
lead_dt <- fread(file.path(RESULTS_DIR, "sle_meta_lead_snps.tsv"))
fg_dt   <- fread(file.path(RESULTS_DIR, "finngen_replication_v2_results_full.tsv"))

cat("Discovery lead variants (after conservative pruning):", nrow(lead_dt), "\n")
cat("FinnGen full results rows:", nrow(fg_dt), "\n\n")

# ── Classify every lead SNP into exactly one replication category ─────────────
# Categories: palindromic_excluded | replicated | not_replicated
fg_status <- fg_dt[, .(SNP, BETA_fg, SE_fg, P_fg, align_status, dir_consistent, replicated)]

# Attach status to all lead SNPs (keep ALL 13)
combined_base <- merge(
  lead_dt[, .(SNP, CHR, BP, A1, A2, BETA_meta=BETA, SE_meta=SE, P_meta=P)],
  fg_status,
  by="SNP", all.x=TRUE
)

# Explicit replication_status column — no ambiguous NA
combined_base[, replication_status := fcase(
  align_status == "palindromic_excluded", "palindromic_excluded",
  replicated == TRUE,                     "replicated",
  replicated == FALSE & dir_consistent == TRUE,  "not_replicated_direction_consistent",
  replicated == FALSE & dir_consistent == FALSE, "not_replicated_discordant_direction",
  is.na(align_status),                    "not_found_in_FinnGen",
  default = "not_replicated"
)]

# ── Attach annotation ─────────────────────────────────────────────────────────
loci_dt <- fread(file.path(RESULTS_DIR, "top_loci_summary_table.tsv"),
                 select = c("SNP","Gene","Region","Known_SLE","Novel"))
combined <- merge(combined_base, loci_dt, by="SNP", all.x=TRUE)

# ── Save full table (all 13 SNPs, all status columns) ────────────────────────
fwrite(combined[order(P_meta)],
       file.path(RESULTS_DIR, "european_replication_concordance_table.tsv"), sep="\t")
cat("Saved: results/european_replication_concordance_table.tsv\n\n")

# ── Print publication-style summary (all 13 rows, status explicit) ────────────
cat("=== European Replication Concordance (publication table) ===\n")
pub_tbl <- combined[order(P_meta), .(
  SNP,
  Locus              = Gene,
  Meta_P             = formatC(P_meta, format="e", digits=2),
  FG_BETA            = round(BETA_fg, 3),
  FG_P               = ifelse(is.na(P_fg), "—", formatC(P_fg, format="e", digits=2)),
  Replication_status = replication_status
)]
print(pub_tbl)

# ── Count breakdown (explicit, no ambiguity) ──────────────────────────────────
n_total       <- nrow(combined)
n_palindromic <- sum(combined$replication_status == "palindromic_excluded")
n_testable    <- n_total - n_palindromic
n_replicated  <- sum(combined$replication_status == "replicated")
n_dir_cons    <- sum(combined$replication_status == "not_replicated_direction_consistent")
n_discordant  <- sum(combined$replication_status == "not_replicated_discordant_direction")
n_not_found   <- sum(combined$replication_status == "not_found_in_FinnGen")

cat(sprintf("\n=== Replication Count Summary ===\n"))
cat(sprintf("  Total lead variants (conservative pruning) : %d\n", n_total))
cat(sprintf("  Excluded — palindromic (strand ambiguous)  : %d  [rs7582694 C/G]\n", n_palindromic))
cat(sprintf("  Testable in FinnGen R12                    : %d\n", n_testable))
cat(sprintf("  Replicated (P < %.4f & concordant dir)    : %d\n",
            0.05 / n_testable, n_replicated))
cat(sprintf("  Not replicated — concordant direction      : %d\n", n_dir_cons))
cat(sprintf("  Not replicated — discordant direction      : %d\n", n_discordant))
cat(sprintf("  Not found in FinnGen                       : %d\n", n_not_found))
cat(sprintf("\n  Headline: %d / %d testable loci replicated in FinnGen R12\n",
            n_replicated, n_testable))

cat("\n=== Step 4 complete ===\n")
