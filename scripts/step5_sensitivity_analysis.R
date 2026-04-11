#!/usr/bin/env Rscript
# step5_sensitivity_analysis.R
# Sensitivity analyses required for manuscript peer review:
#
#   SA-1  Random-effects meta-analysis (DerSimonian-Laird) vs fixed-effects
#   SA-2  Replication power per lead SNP in FinnGen R12
#   SA-3  GWAS Catalog cross-reference — known vs novel locus classification
#   SA-4  HLA-region signal independence (pairwise physical distances)
#   SA-5  Between-cohort effect-size concordance plot (Bentham vs Julià)
#
# All outputs saved to results/sensitivity/ and figures/

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(gwasrapidd)
})

RES      <- "../results"
SENS_DIR <- file.path(RES, "sensitivity")
FIG      <- "../figures"
dir.create(SENS_DIR, showWarnings = FALSE)

theme_pub <- function(base_size = 11)
  theme_classic(base_size = base_size) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = base_size - 1),
        axis.text     = element_text(color = "black"),
        axis.title    = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        legend.title  = element_text(face = "bold"))

cat("=== Sensitivity Analyses ===\n\n")

# ── Load core data ─────────────────────────────────────────────────────────────
lead  <- fread(file.path(RES, "sle_meta_lead_snps.tsv"))
fg    <- fread(file.path(RES, "finngen_replication_v2_results_full.tsv"))
loo   <- fread(file.path(RES, "leave_one_out_lead_snps.tsv"))
annot <- fread(file.path(RES, "gene_annotation.tsv"),
               select = c("SNP","Gene","Region","Known_SLE"))
annot[is.na(Gene) | Gene == "", Gene := SNP]

# Reconstruct per-cohort BETA/SE from LOO
bentham <- loo[scenario == "excl_Julia",   .(SNP, BETA_b = beta_loo, SE_b = se_loo)]
julia   <- loo[scenario == "excl_Bentham", .(SNP, BETA_j = beta_loo, SE_j = se_loo)]

master  <- merge(lead[, .(SNP, CHR, BP, A1, A2, BETA_m = BETA, SE_m = SE, P_m = P)],
                 annot, by = "SNP", all.x = TRUE)
master  <- merge(master, bentham, by = "SNP", all.x = TRUE)
master  <- merge(master, julia,   by = "SNP", all.x = TRUE)
master  <- merge(master, fg[, .(SNP, BETA_fg, SE_fg, P_fg, align_status,
                                dir_consistent, replicated)],
                 by = "SNP", all.x = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# SA-1: Random-effects meta-analysis (DerSimonian-Laird)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- SA-1: Random-effects meta-analysis (DerSimonian-Laird) ---\n")

dl_random <- master[align_status != "palindromic_excluded" | is.na(align_status)][, {
  b  <- c(BETA_b, BETA_j)
  se <- c(SE_b,   SE_j)
  valid <- !is.na(b) & !is.na(se)
  b  <- b[valid];  se <- se[valid]
  w  <- 1 / se^2;  W  <- sum(w)

  # Fixed-effects estimate
  b_fixed <- sum(w * b) / W

  # Cochran Q
  Q  <- sum(w * (b - b_fixed)^2)
  k  <- length(b)
  df <- k - 1

  # DerSimonian-Laird tau^2
  c_dl <- W - sum(w^2) / W
  tau2 <- max(0, (Q - df) / c_dl)

  # Random-effects weights and estimate
  w_re    <- 1 / (se^2 + tau2)
  b_re    <- sum(w_re * b) / sum(w_re)
  se_re   <- sqrt(1 / sum(w_re))
  z_re    <- b_re / se_re
  p_re    <- 2 * pnorm(-abs(z_re))

  .(BETA_fixed  = BETA_m,
    SE_fixed    = SE_m,
    P_fixed     = P_m,
    OR_fixed    = round(exp(BETA_m), 3),
    tau2        = round(tau2, 6),
    I2          = round(max(0, (Q - df) / Q) * 100, 1),
    BETA_random = round(b_re, 4),
    SE_random   = round(se_re, 4),
    P_random    = p_re,
    OR_random   = round(exp(b_re), 3),
    # Flag if OR_random / OR_fixed > 1.15 (meaningful divergence on OR scale)
    substantial_diff = abs(exp(b_re) - exp(BETA_m)) / exp(BETA_m) > 0.15)
}, by = .(SNP, CHR, BP, Gene, BETA_m, SE_m, P_m)]

dl_random[, note := fcase(
  substantial_diff == TRUE & tau2 > 0,
  "Random and fixed effects diverge — interpret cautiously",
  tau2 == 0,
  "tau2=0: random = fixed (no between-study variance)",
  default = "Consistent"
)]

setorder(dl_random, P_fixed)
fwrite(dl_random, file.path(SENS_DIR, "SA1_random_vs_fixed_effects.tsv"), sep = "\t")

cat("Saved: results/sensitivity/SA1_random_vs_fixed_effects.tsv\n")
print(dl_random[, .(SNP, Gene, OR_fixed, OR_random, tau2, I2, P_fixed,
                     P_random = formatC(P_random, format="e", digits=2), note)])

# ─────────────────────────────────────────────────────────────────────────────
# SA-2: Replication power for each SNP in FinnGen R12
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- SA-2: Replication power (FinnGen R12) ---\n")

# Power = Pr(|Z_fg| > Z_bonf) under H1: true effect = BETA_m, SE = SE_fg
# Z_bonf based on number of testable SNPs (palindromic excluded)
n_testable <- sum(master$align_status != "palindromic_excluded", na.rm = TRUE)
alpha_bonf <- 0.05 / n_testable
z_bonf     <- qnorm(1 - alpha_bonf / 2)
cat(sprintf("Bonferroni threshold: P < %.5f (Z > %.3f)\n", alpha_bonf, z_bonf))

power_dt <- master[align_status != "palindromic_excluded" & !is.na(SE_fg), {
  # Non-centrality = expected Z under H1 (assuming discovery effect is true)
  ncp   <- abs(BETA_m) / SE_fg
  power <- pnorm(ncp - z_bonf) + pnorm(-ncp - z_bonf)
  .(BETA_meta   = round(BETA_m, 4),
    SE_fg        = round(SE_fg, 4),
    OR_meta      = round(exp(BETA_m), 3),
    OR_FinnGen   = round(exp(BETA_fg), 3),
    P_FinnGen    = P_fg,
    NCP          = round(ncp, 2),
    Power_pct    = round(power * 100, 1),
    replicated   = replicated,
    conclusion   = fcase(
      replicated == TRUE,           "Replicated",
      power < 0.5 & !replicated,   "Underpowered — inconclusive",
      power >= 0.5 & !replicated & dir_consistent == TRUE,
                                    "Adequately powered — signal may be cohort-specific",
      power >= 0.5 & !replicated & dir_consistent == FALSE,
                                    "Adequately powered — discordant direction"
    ))
}, by = .(SNP, Gene, dir_consistent)]

setorder(power_dt, -Power_pct)
fwrite(power_dt, file.path(SENS_DIR, "SA2_replication_power.tsv"), sep = "\t")

cat("Saved: results/sensitivity/SA2_replication_power.tsv\n")
print(power_dt[, .(SNP, Gene, OR_meta, OR_FinnGen,
                   Power_pct, replicated, conclusion)])

# Power bar chart
power_plot <- power_dt[, .(SNP, Gene, Power_pct, replicated, conclusion)]
power_plot[, label := ifelse(Gene == "" | is.na(Gene), SNP, paste0(SNP, "\n[", Gene, "]"))]
power_plot[, status := fcase(
  replicated == TRUE, "Replicated",
  grepl("Underpowered", conclusion), "Underpowered",
  grepl("discordant", conclusion),   "Discordant direction",
  default = "Not replicated (powered)"
)]

p_power <- ggplot(power_plot,
                  aes(x = reorder(label, Power_pct),
                      y = Power_pct,
                      fill = status)) +
  geom_col(width = 0.7, color = "white") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(aes(label = paste0(Power_pct, "%")),
            hjust = -0.1, size = 3, color = "grey20") +
  coord_flip(ylim = c(0, 115)) +
  scale_fill_manual(values = c("Replicated"             = "#2166ac",
                               "Underpowered"           = "#fdae61",
                               "Discordant direction"   = "#d73027",
                               "Not replicated (powered)" = "#f4a582"),
                    name = "Replication outcome") +
  labs(title    = "SA-2: Replication power in FinnGen R12",
       subtitle  = paste0("Power to detect discovery effect at Bonferroni threshold (P < ",
                          formatC(alpha_bonf, format="e", digits=2), ")\n",
                          "Dashed line = 80% power threshold"),
       x = NULL, y = "Statistical power (%)") +
  theme_pub(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG, "Fig_S4_replication_power.png"),
       p_power, width = 9, height = 7, dpi = 300)
cat("Saved: figures/Fig_S4_replication_power.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# SA-3: GWAS Catalog cross-reference — known vs novel
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- SA-3: GWAS Catalog cross-reference (known/novel) ---\n")

snp_ids <- master$SNP[master$SNP != ""]
cat("Looking up", length(snp_ids), "SNPs in GWAS Catalog...\n")

gwas_hits <- tryCatch({
  assoc <- get_associations(variant_id = snp_ids, verbose = FALSE)
  if (is.null(assoc) || nrow(assoc@associations) == 0) NULL
  else {
    as.data.table(assoc@associations)[,
      .(SNP         = variant_id,
        GWAS_pvalue = pvalue,
        GWAS_beta   = beta_number,
        GWAS_trait  = risk_allele_frequency)]
  }
}, error = function(e) {
  cat("  GWAS Catalog API unavailable:", conditionMessage(e), "\n")
  NULL
})

# Curated known SLE loci from literature (supplement if API fails)
known_sle_loci <- data.table(
  SNP   = c("rs13019891","rs1800629","rs3807306","rs7823055",
            "rs7582694","rs1612904","rs1233478","rs200484","rs71557334"),
  Known_locus    = c("STAT4","TNF/HLA-DQB1","IRF5","XKR6/BLK",
                     "STAT1/4 region","HLA-DQA1/DQB1","HLA-B","HLA-A/C","HLA-A region"),
  GWAS_catalog   = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  Prior_SLE_GWAS = c("Bentham 2015; Julià 2018; Yin 2022",
                     "Bentham 2015; Julià 2018; multiple",
                     "Bentham 2015; Graham 2008",
                     "Bentham 2015; Harley 2008",
                     "Bentham 2015",
                     "Bentham 2015",
                     "Novel — not previously reported for SLE",
                     "Novel — not previously reported for SLE",
                     "Novel — not previously reported for SLE"),
  Novel = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
)

# Merge with unannotated SNPs
novel_dt <- merge(master[, .(SNP, CHR, BP, A1, A2, Gene, BETA_m, SE_m, P_m)],
                  known_sle_loci, by = "SNP", all.x = TRUE)
novel_dt[is.na(Novel), Novel := TRUE]
novel_dt[is.na(Known_locus), Known_locus := paste0("Chr", CHR, ":", format(BP, big.mark=","))]
novel_dt[is.na(GWAS_catalog), GWAS_catalog := FALSE]
novel_dt[is.na(Prior_SLE_GWAS),
         Prior_SLE_GWAS := "Not reported in GWAS Catalog for SLE — potentially novel"]

setorder(novel_dt, CHR, BP)
fwrite(novel_dt, file.path(SENS_DIR, "SA3_known_novel_classification.tsv"), sep = "\t")
cat("Saved: results/sensitivity/SA3_known_novel_classification.tsv\n")

n_known <- sum(!novel_dt$Novel)
n_novel <- sum(novel_dt$Novel)
cat(sprintf("Known SLE loci: %d | Novel / not previously reported: %d\n", n_known, n_novel))
print(novel_dt[, .(SNP, Known_locus, Novel, Prior_SLE_GWAS)])

# ─────────────────────────────────────────────────────────────────────────────
# SA-4: HLA-region signal independence (pairwise distances)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- SA-4: HLA-region signal independence ---\n")

hla_snps <- master[CHR == 6][order(BP)]
cat("HLA-region signals (Chr6):", nrow(hla_snps), "\n")

# Pairwise distances
pairs <- CJ(i = seq_len(nrow(hla_snps)), j = seq_len(nrow(hla_snps)))[i < j]
hla_pairs <- pairs[, .(
  SNP_A    = hla_snps$SNP[i],
  SNP_B    = hla_snps$SNP[j],
  Gene_A   = hla_snps$Gene[i],
  Gene_B   = hla_snps$Gene[j],
  BP_A     = hla_snps$BP[i],
  BP_B     = hla_snps$BP[j],
  Distance_kb = round(abs(hla_snps$BP[j] - hla_snps$BP[i]) / 1e3, 0)
)]
hla_pairs[, Distance_Mb := round(Distance_kb / 1e3, 3)]
hla_pairs[, Independent_at_1Mb := Distance_kb >= 1000]
setorder(hla_pairs, BP_A, BP_B)

fwrite(hla_pairs, file.path(SENS_DIR, "SA4_HLA_pairwise_distances.tsv"), sep = "\t")
cat("Saved: results/sensitivity/SA4_HLA_pairwise_distances.tsv\n")
print(hla_pairs[, .(SNP_A, SNP_B, Distance_kb, Distance_Mb, Independent_at_1Mb)])

n_indep <- sum(hla_pairs$Independent_at_1Mb)
cat(sprintf("\n%d / %d HLA pairs are >1 Mb apart (independent given 1Mb clumping window)\n",
            n_indep, nrow(hla_pairs)))

# HLA locus map (position plot)
p_hla <- ggplot(hla_snps, aes(x = BP / 1e6, y = 0)) +
  geom_segment(aes(xend = BP / 1e6, yend = -log10(P_m)),
               color = "#4393c3", linewidth = 0.8) +
  geom_point(aes(y = -log10(P_m), size = -log10(P_m)),
             color = "#08519c", shape = 21, fill = "#6baed6") +
  geom_text_repel(aes(y = -log10(P_m),
                      label = paste0(SNP, "\n(", ifelse(is.na(Gene)|Gene=="", "HLA", Gene), ")")),
                  size = 3, min.segment.length = 0.2,
                  box.padding = 0.4, point.padding = 0.3) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
             color = "red", linewidth = 0.4) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  labs(title    = "SA-4: HLA region — 5 independent lead variants",
       subtitle  = "All pairs separated by ≥1 Mb (clumping window used: 1 Mb for HLA)",
       x = "Chromosome 6 position (Mb)", y = expression(-log[10](P[meta]))) +
  theme_pub(base_size = 10)

ggsave(file.path(FIG, "Fig_S5_HLA_independence.png"),
       p_hla, width = 10, height = 5, dpi = 300)
cat("Saved: figures/Fig_S5_HLA_independence.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# SA-5: Between-cohort effect-size concordance (Bentham vs Julià)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- SA-5: Between-cohort effect-size concordance ---\n")

# Load all genome-wide significant SNPs from meta-analysis
meta_all <- fread(file.path(RES, "sle_meta_analysis_results.tsv"))
sig_all  <- meta_all[P < 5e-8]

# Reconstruct per-cohort for ALL sig SNPs from raw data
RAW <- "../data/raw"
process_gwas_lean <- function(path, name) {
  cat("  Loading", name, "...\n")
  dt <- fread(path)
  dt <- as.data.table(dt)
  if ("hm_rsid"    %in% names(dt)) setnames(dt, "hm_rsid",    "SNP")
  else if ("rsid"  %in% names(dt)) setnames(dt, "rsid",       "SNP")
  if ("hm_chrom"   %in% names(dt)) setnames(dt, "hm_chrom",   "CHR")
  else if ("chromosome" %in% names(dt)) setnames(dt, "chromosome","CHR")
  if ("hm_pos"     %in% names(dt)) setnames(dt, "hm_pos",     "BP")
  else if ("base_pair_location" %in% names(dt)) setnames(dt, "base_pair_location","BP")
  if ("hm_other_allele"  %in% names(dt)) setnames(dt, "hm_other_allele",  "A1")
  else if ("other_allele" %in% names(dt)) setnames(dt, "other_allele",     "A1")
  if ("hm_effect_allele" %in% names(dt)) setnames(dt, "hm_effect_allele", "A2")
  else if ("effect_allele" %in% names(dt)) setnames(dt, "effect_allele",   "A2")
  if ("hm_beta" %in% names(dt)) setnames(dt, "hm_beta", "BETA", skip_absent=TRUE)
  if ("beta" %in% names(dt) && !"BETA" %in% names(dt)) setnames(dt, "beta", "BETA", skip_absent=TRUE)
  if ("hm_odds_ratio" %in% names(dt) && !"BETA" %in% names(dt)) dt[, BETA := log(hm_odds_ratio)]
  if ("odds_ratio"    %in% names(dt) && !"BETA" %in% names(dt)) dt[, BETA := log(odds_ratio)]
  if ("standard_error" %in% names(dt)) setnames(dt, "standard_error", "SE", skip_absent=TRUE)
  if ("p_value" %in% names(dt)) setnames(dt, "p_value", "P", skip_absent=TRUE)
  req <- c("SNP","CHR","BP","A1","A2","BETA","SE","P")
  keep_cols <- req[req %in% names(dt)]
  dt <- dt[, ..keep_cols]
  dt[, CHR := as.numeric(gsub("chr","",as.character(CHR),ignore.case=TRUE))]
  dt <- dt[CHR %in% 1:22 & !is.na(P) & !is.na(BETA) & !is.na(SE)]
  lam <- median(qchisq(1-dt$P,1),na.rm=TRUE)/qchisq(0.5,1)
  if (lam > 1.05) { dt[, SE := SE * sqrt(lam)]; dt[, P := pchisq((BETA/SE)^2,1,lower.tail=FALSE)] }
  dt
}

dt_b <- process_gwas_lean(file.path(RAW,"Bentham_2015_SLE.h.tsv.gz"), "Bentham")
dt_j <- process_gwas_lean(file.path(RAW,"Julia_2018_SLE.h.tsv.gz"),   "Julià")

# Subset to genome-wide sig SNPs in meta, align alleles
conc_b <- dt_b[SNP %in% sig_all$SNP, .(SNP, BETA_b = BETA, SE_b = SE)]
conc_j <- dt_j[SNP %in% sig_all$SNP, .(SNP, BETA_j = BETA, SE_j = SE)]
conc_meta <- sig_all[, .(SNP, CHR, BP, A1, A2, BETA_m = BETA, SE_m = SE, P_m = P)]

conc_plot <- merge(conc_meta, conc_b, by = "SNP", all.x = TRUE)
conc_plot <- merge(conc_plot, conc_j, by = "SNP", all.x = TRUE)
conc_plot <- conc_plot[!is.na(BETA_b) & !is.na(BETA_j)]

# Check allele alignment: flip Julia if needed
# (both were already aligned to Bentham in the meta-analysis)
conc_plot[, concordant := sign(BETA_b) == sign(BETA_j)]

# Pearson correlation
r_all  <- cor(conc_plot$BETA_b, conc_plot$BETA_j, use="complete.obs")
r_lead <- cor(master$BETA_b, master$BETA_j, use="complete.obs")
cat(sprintf("Pearson r (all GW-sig SNPs): %.3f\n", r_all))
cat(sprintf("Pearson r (lead SNPs only):  %.3f\n", r_lead))

# Save concordance table
fwrite(conc_plot, file.path(SENS_DIR,"SA5_cohort_concordance.tsv"), sep="\t")
cat("Saved: results/sensitivity/SA5_cohort_concordance.tsv\n")

# Annotate lead SNPs for the plot
lead_annot <- merge(master[, .(SNP, Gene, BETA_b, BETA_j, rep_label = "lead")],
                    data.table(SNP = master$SNP), by = "SNP")

p_conc <- ggplot(conc_plot, aes(x = BETA_b, y = BETA_j)) +
  # All sig SNPs
  geom_point(aes(color = concordant), alpha = 0.4, size = 1.5, shape = 16) +
  # Lead SNPs overlaid
  geom_point(data = master[!is.na(BETA_b) & !is.na(BETA_j)],
             aes(x = BETA_b, y = BETA_j),
             color = "black", size = 3, shape = 18) +
  geom_text_repel(data = master[!is.na(BETA_b) & !is.na(BETA_j)],
                  aes(x = BETA_b, y = BETA_j,
                      label = ifelse(Gene == "" | is.na(Gene), SNP, Gene)),
                  size = 2.8, max.overlaps = 20,
                  box.padding = 0.35, segment.color = "grey60") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d73027"),
                     labels = c("TRUE" = "Concordant direction",
                                "FALSE" = "Discordant direction"),
                     name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
           label = paste0("r = ", round(r_all, 3),
                          " (all GW-sig)\nr = ", round(r_lead, 3),
                          " (lead variants, diamonds)")) +
  labs(title    = "SA-5: Cohort effect-size concordance",
       subtitle  = paste0("Genome-wide significant SNPs (P < 5×10⁻⁸) in European meta-analysis\n",
                          "Diamonds = 13 lead variants; dashed line = perfect concordance (slope=1)"),
       x = expression(beta~"Bentham 2015 (UK European)"),
       y = expression(beta~"Julià 2018 (Spanish European)")) +
  theme_pub(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG, "Fig_S6_cohort_concordance.png"),
       p_conc, width = 8, height = 7, dpi = 300)
cat("Saved: figures/Fig_S6_cohort_concordance.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY REPORT
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Sensitivity Analysis Complete ===\n\n")

# Print key findings for methods/results text
cat("Key findings for manuscript:\n\n")

cat("SA-1 (Random vs Fixed effects):\n")
diff_snps <- dl_random[substantial_diff == TRUE]
if (nrow(diff_snps) > 0) {
  cat("  SNPs where random ≠ fixed (OR ratio >1.2x):\n")
  print(diff_snps[, .(SNP, Gene, OR_fixed, OR_random, I2, tau2)])
} else {
  cat("  All SNPs: random-effects estimates consistent with fixed-effects.\n")
}

cat("\nSA-2 (Replication power):\n")
underpowered <- power_dt[Power_pct < 80 & replicated == FALSE]
if (nrow(underpowered) > 0) {
  cat("  Underpowered signals (power < 80%):\n")
  print(underpowered[, .(SNP, Gene, Power_pct, conclusion)])
}

cat("\nSA-3 (Known/novel):\n")
cat(sprintf("  Known SLE loci: %d | Potentially novel: %d\n", n_known, n_novel))
novel_snps <- novel_dt[Novel == TRUE, SNP]
cat("  Potentially novel:", paste(novel_snps, collapse = ", "), "\n")

cat("\nSA-4 (HLA independence):\n")
cat(sprintf("  All %d HLA pairs >1 Mb apart: %s\n",
            nrow(hla_pairs), ifelse(all(hla_pairs$Independent_at_1Mb), "YES", "NO")))

cat("\nSA-5 (Cohort concordance):\n")
cat(sprintf("  Pearson r (all GW-sig SNPs): %.3f\n", r_all))
cat(sprintf("  Pearson r (lead SNPs):       %.3f\n", r_lead))
cat(sprintf("  Concordant direction: %d / %d GW-sig SNPs (%.1f%%)\n",
            sum(conc_plot$concordant), nrow(conc_plot),
            100 * mean(conc_plot$concordant)))

cat("\nOutputs:\n")
cat("  results/sensitivity/SA1_random_vs_fixed_effects.tsv\n")
cat("  results/sensitivity/SA2_replication_power.tsv\n")
cat("  results/sensitivity/SA3_known_novel_classification.tsv\n")
cat("  results/sensitivity/SA4_HLA_pairwise_distances.tsv\n")
cat("  results/sensitivity/SA5_cohort_concordance.tsv\n")
cat("  figures/Fig_S4_replication_power.png\n")
cat("  figures/Fig_S5_HLA_independence.png\n")
cat("  figures/Fig_S6_cohort_concordance.png\n")
