library(data.table)
library(dplyr)
library(ggplot2)

RAW_DIR     <- "../data/raw"
RESULTS_DIR <- "../results"
FIG_DIR     <- "../figures"
dir.create(FIG_DIR, showWarnings = FALSE)

cat("=== STEP 1: Heterogeneity, Forest Plots, Per-Cohort Consistency, Leave-One-Out (European-only) ===\n\n")

# ── Load processed datasets by re-running the munging logic from the main script ──
process_gwas <- function(file_path, name) {
  cat("Loading", name, "...\n")
  dt <- fread(file_path)
  
  if ("hm_rsid" %in% names(dt)) {
    setnames(dt, old = "hm_rsid", new = "SNP")
  } else if ("rsid" %in% names(dt)) {
    setnames(dt, old = "rsid", new = "SNP")
  } else if ("variant_id" %in% names(dt)) {
    setnames(dt, old = "variant_id", new = "SNP")
  }
  
  if("hm_chrom" %in% names(dt)){ setnames(dt, old="hm_chrom", new="CHR") }
  else if("chromosome" %in% names(dt)){ setnames(dt, old="chromosome", new="CHR") }
  
  if("hm_pos" %in% names(dt)){ setnames(dt, old="hm_pos", new="BP") }
  else if("base_pair_location" %in% names(dt)){ setnames(dt, old="base_pair_location", new="BP") }
  
  if ("hm_other_allele" %in% names(dt)) setnames(dt, "hm_other_allele", "A1")
  else if ("other_allele" %in% names(dt)) setnames(dt, "other_allele", "A1")
  else if ("AlleleB" %in% names(dt)) setnames(dt, "AlleleB", "A1")
  
  if ("hm_effect_allele" %in% names(dt)) setnames(dt, "hm_effect_allele", "A2")
  else if ("effect_allele" %in% names(dt)) setnames(dt, "effect_allele", "A2")
  else if ("A1leleA" %in% names(dt)) setnames(dt, "A1leleA", "A2")
  
  if ("hm_beta" %in% names(dt)) setnames(dt, "hm_beta", "BETA", skip_absent=TRUE)
  if ("beta" %in% names(dt) && !"BETA" %in% names(dt)) setnames(dt, "beta", "BETA", skip_absent=TRUE)
  if ("hm_odds_ratio" %in% names(dt) && !"BETA" %in% names(dt)) dt[, BETA := log(hm_odds_ratio)]
  if ("odds_ratio" %in% names(dt) && !"BETA" %in% names(dt)) dt[, BETA := log(odds_ratio)]
  if ("OR" %in% names(dt) && !"BETA" %in% names(dt)) dt[, BETA := log(OR)]
  
  if ("standard_error" %in% names(dt)) setnames(dt, "standard_error", "SE", skip_absent=TRUE)
  if ("p_value" %in% names(dt)) setnames(dt, "p_value", "P", skip_absent=TRUE)
  if ("P" %in% names(dt) & !"P" %in% names(dt)) setnames(dt, "P", "P", skip_absent=TRUE)

  req_cols <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
  present_cols <- intersect(req_cols, names(dt))
  missing_cols <- setdiff(req_cols, present_cols)
  
  if(length(missing_cols) > 0){
     if ("SE" %in% missing_cols && "P" %in% names(dt) && "BETA" %in% names(dt)) {
         dt[, SE := abs(BETA) / qnorm(1 - P/2)]
         present_cols <- c(present_cols, "SE")
     }
  }

  dt <- dt[, ..present_cols]
  dt[, CHR := as.numeric(gsub("chr","",as.character(CHR),ignore.case=TRUE))]
  dt <- dt[CHR %in% 1:22]
  dt[, A1 := toupper(A1)]; dt[, A2 := toupper(A2)]
  dt <- dt[!is.na(P) & !is.na(SE) & !is.na(BETA)]
  
  # GC correction
  lambda <- median(qchisq(1-dt$P,1),na.rm=TRUE)/qchisq(0.5,1)
  cat("  Lambda:", round(lambda,3), "\n")
  if (lambda > 1.05) {
    dt[, SE := SE * sqrt(lambda)]
    dt[, P  := pchisq((BETA/SE)^2, df=1, lower.tail=FALSE)]
  }
  dt
}

dt_bentham <- process_gwas(file.path(RAW_DIR,"Bentham_2015_SLE.h.tsv.gz"), "Bentham")
dt_julia   <- process_gwas(file.path(RAW_DIR,"Julia_2018_Spain_Only.txt"), "Julià")
# Yin 2022 (East Asian) excluded — European-ancestry analysis only

# Load lead SNPs
lead_dt <- fread(file.path(RESULTS_DIR, "sle_meta_lead_snps.tsv"))
lead_snps <- lead_dt$SNP
cat("\nLead SNPs:", paste(lead_snps, collapse=", "), "\n")

# ── Subset each cohort to lead SNPs ──
b_lead <- dt_bentham[SNP %in% lead_snps, .(SNP, BETA_b=BETA, SE_b=SE, P_b=P)]
j_lead <- dt_julia  [SNP %in% lead_snps, .(SNP, BETA_j=BETA, SE_j=SE, P_j=P)]
meta_lead <- lead_dt[, .(SNP, CHR, BP, A1, A2, BETA_m=BETA, SE_m=SE, P_m=P)]

all_cohorts <- Reduce(function(a,b) merge(a,b,by="SNP",all=TRUE),
                      list(meta_lead, b_lead, j_lead))

# ── 1a. Per-cohort direction of effect ──
cat("\n--- 1a. Per-cohort direction of effect ---\n")
pcoh <- all_cohorts[, .(SNP, CHR, BP,
                        Dir_Bentham = sign(BETA_b),
                        Dir_Julia   = sign(BETA_j),
                        Dir_Meta    = sign(BETA_m))]
pcoh[, Consistent := (Dir_Bentham == Dir_Julia)]
fwrite(pcoh, file.path(RESULTS_DIR,"per_cohort_effects.tsv"), sep="\t")
print(pcoh)

# ── 1b. Cochran's Q and I² for lead SNPs ──
cat("\n--- 1b. Heterogeneity (Cochran Q, I²) for lead SNPs ---\n")
het <- all_cohorts %>%
  rowwise() %>%
  mutate(
    betas   = list(c(BETA_b, BETA_j)),
    ses     = list(c(SE_b,   SE_j)),
    weights = list(1/c(SE_b, SE_j)^2),
    W_total = sum(unlist(weights), na.rm=TRUE),
    meta_b  = sum(unlist(weights)*unlist(betas), na.rm=TRUE) / W_total,
    Q       = sum(unlist(weights)*(unlist(betas)-meta_b)^2, na.rm=TRUE),
    df      = sum(!is.na(unlist(betas))) - 1,
    p_Q     = pchisq(Q, df=df, lower.tail=FALSE),
    I2      = round(max(0, (Q-df)/Q)*100, 1)
  ) %>%
  select(SNP, Q, df, p_Q, I2)

fwrite(as.data.table(het), file.path(RESULTS_DIR,"heterogeneity_lead_snps.tsv"), sep="\t")
print(het)

# ── Heterogeneity interpretation note ────────────────────────────────────────
# With only 2 discovery cohorts Cochran's Q has df = 1. This gives very low
# power to detect heterogeneity at modest I², but makes high I² values (e.g.
# >85%) highly significant statistically. For HLA-region loci specifically,
# high I² is biologically plausible: the HLA harbours extensive long-range LD
# and haplotype frequency differences even between North European (Bentham;
# UK) and South European (Julià; Spanish) populations. Heterogeneity here
# should be interpreted cautiously — it may reflect true locus complexity and
# population-specific haplotype structures rather than analytical artefact.
# Caution is advised before concluding heterogeneity invalidates a signal;
# the replication result in the independent Finnish FinnGen cohort is the
# primary arbiter of signal validity.
het_dt <- as.data.table(het)
hla_snps <- c("rs1800629","rs1612904","rs1233478","rs200484","rs71557334")
hla_het  <- het_dt[SNP %in% hla_snps]
if (nrow(hla_het) > 0) {
  cat("\n--- HLA-region heterogeneity (interpret cautiously; df = 1) ---\n")
  cat("High I² in HLA is plausible given population-specific haplotype\n")
  cat("frequency differences between UK (Bentham) and Spanish (Julià) cohorts.\n")
  cat("This does NOT necessarily indicate analytical artefact.\n\n")
  print(hla_het[, .(SNP, I2, p_Q)])
}

# ── 1c. Forest plots for each lead SNP ──
cat("\n--- 1c. Forest plots ---\n")
for (snp_id in lead_snps) {
  row <- all_cohorts[SNP == snp_id]
  if (nrow(row) == 0) next

  plot_df <- data.frame(
    cohort = c("Bentham","Julià","■ Meta"),
    beta   = c(row$BETA_b, row$BETA_j, row$BETA_m),
    se     = c(row$SE_b,   row$SE_j,   row$SE_m)
  ) %>%
    filter(!is.na(beta)) %>%
    mutate(
      lower = beta - 1.96*se,
      upper = beta + 1.96*se,
      is_meta = cohort == "■ Meta"
    )

  p <- ggplot(plot_df, aes(x=beta, y=cohort, color=is_meta)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey50") +
    geom_point(aes(size=is_meta)) +
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.2) +
    scale_color_manual(values=c("FALSE"="steelblue","TRUE"="darkred"), guide="none") +
    scale_size_manual(values=c("FALSE"=2,"TRUE"=3), guide="none") +
    labs(title=paste0("Forest Plot: ", snp_id),
         x="Effect size (β ± 95% CI)", y=NULL) +
    theme_minimal(base_size=11) +
    theme(plot.title=element_text(face="bold"))

  ggsave(file.path(FIG_DIR, paste0("forest_",snp_id,".png")), p, width=6, height=3.5, dpi=150)
}
cat("Forest plots saved to figures/\n")

# ── 1d. Leave-one-out meta-analysis (for lead SNPs) ──
cat("\n--- 1d. Leave-one-cohort-out ---\n")

ivw_meta <- function(betas, ses) {
  w <- 1/ses^2
  mb <- sum(w*betas, na.rm=TRUE)/sum(w, na.rm=TRUE)
  mse <- sqrt(1/sum(w, na.rm=TRUE))
  mz  <- mb/mse
  mp  <- 2*pnorm(-abs(mz))
  c(beta=mb, se=mse, p=mp)
}

loo_results <- lapply(list(
  list(name="excl_Bentham", b=all_cohorts$BETA_j, se_b=all_cohorts$SE_j),
  list(name="excl_Julia",   b=all_cohorts$BETA_b, se_b=all_cohorts$SE_b)
), function(cfg) {
  dt_loo <- all_cohorts[, .(SNP)]
  # With 2 discovery cohorts, LOO = single remaining cohort effect
  p_val <- 2 * pnorm(-abs(cfg$b / cfg$se_b))
  dt_loo[, `:=`(beta_loo=cfg$b, se_loo=cfg$se_b, p_loo=p_val, scenario=cfg$name)]
  dt_loo
})

loo_dt <- rbindlist(loo_results)
fwrite(loo_dt, file.path(RESULTS_DIR,"leave_one_out_lead_snps.tsv"), sep="\t")
print(loo_dt)
cat("\n=== Step 1 complete ===\n")
