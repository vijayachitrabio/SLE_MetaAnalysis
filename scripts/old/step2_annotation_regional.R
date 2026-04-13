library(data.table)
library(dplyr)
library(ggplot2)

RESULTS_DIR <- "../results"
FIG_DIR     <- "../figures"
dir.create(FIG_DIR, showWarnings = FALSE)

cat("=== STEP 2: Gene Annotation & Regional Association Plots ===\n\n")

lead_dt <- fread(file.path(RESULTS_DIR, "sle_meta_lead_snps.tsv"))
meta_dt <- fread(file.path(RESULTS_DIR, "sle_meta_analysis_results.tsv"))

# ── 2a: Gene annotation (nearest gene lookup) ──
# Known SLE gene map (curated list from literature)
known_sle_genes <- c(
  "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-A","HLA-B","HLA-C",
  "HLA-DRB5","HLA-DPB1","STAT4","IRF5","BLK","ITGAM","TNFSF4",
  "PTPN22","BANK1","IRAK1","MECP2","TNFAIP3","IKBKE","CFB","C4A","C4B",
  "TREX1","DNASE1L3","XKR6","ATG5","RASGRP3","UHRF1BP1","IKZF2","ETS1"
)

# Nearest gene per lead SNP (based on known SLE loci / HLA region coordinates)
gene_map <- data.table(
  SNP       = c("rs1800629","rs9268671","rs13019891","rs2647046","rs3131783",
                "rs7574865","rs652888","rs2524129","rs3807306","rs2517600",
                "rs11961853","rs204989","rs7823055","rs71557334","rs9257803","rs200484"),
  Gene      = c("TNF/HLA","HLA-DQA1","STAT4","HLA-DQB1","MICB/HLA",
                "STAT4","HLA-C","HCP5/HLA","IRF5","C4B/HLA",
                "HLA-A","HLA-DMA","XKR6/BLK","HLA-upstream","HLA-upstream","HLA-upstream"),
  Region    = c("6p21.3","6p21.3","2q32","6p21.3","6p21.3",
                "2q32","6p21.3","6p21.3","7q32","6p21.3",
                "6p22.1","6p21.3","8p23","6p22.2","6p22.3","6p22.3"),
  Known_SLE = c(TRUE,TRUE,TRUE,TRUE,TRUE,
                TRUE,TRUE,TRUE,TRUE,TRUE,
                TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)
)

gene_annot <- merge(lead_dt, gene_map, by="SNP", all.x=TRUE)
fwrite(gene_annot, file.path(RESULTS_DIR,"gene_annotation.tsv"), sep="\t")
cat("Gene annotation:\n")
print(gene_annot[, .(SNP, CHR, BP, Gene, Region, Known_SLE, BETA, SE, P)])

# ── 2b: Functional annotation (variant class from allele change) ──
cat("\n--- 2b. Functional annotation (variant class) ---\n")
func_annot <- gene_annot %>%
  mutate(
    allele_change = paste0(A1,">",A2),
    variant_class = case_when(
      nchar(A1)==1 & nchar(A2)==1 ~ "SNV",
      nchar(A1) > nchar(A2)      ~ "Deletion",
      nchar(A1) < nchar(A2)      ~ "Insertion",
      TRUE ~ "Other"
    ),
    # Known SLE eQTL genes from immune cells (curated from eQTLGen / GTEx)
    eQTL_gene = case_when(
      SNP == "rs1800629"  ~ "TNF (whole blood, immune cells)",
      SNP == "rs9268671"  ~ "HLA-DQA1 (whole blood)",
      SNP == "rs2647046"  ~ "HLA-DQB1 (whole blood)",
      SNP == "rs7574865"  ~ "STAT4 (whole blood, T cells)",
      SNP == "rs3807306"  ~ "IRF5 (monocytes, pDCs)",
      SNP == "rs13019891" ~ "STAT4 (T cells)",
      TRUE ~ "Not curated"
    )
  )

fwrite(as.data.table(func_annot), file.path(RESULTS_DIR,"functional_annotation.tsv"), sep="\t")
cat("Functional annotation saved.\n")
print(func_annot[, c("SNP","Gene","allele_change","variant_class","eQTL_gene")])

# ── 2c: Regional association plots (LocusZoom-style) ──
cat("\n--- 2c. Regional association plots ---\n")

plot_region <- function(chr, center_bp, label, window_mb=1) {
  window_bp <- window_mb * 1e6
  region_dt <- meta_dt[CHR==chr & BP >= (center_bp - window_bp) & BP <= (center_bp + window_bp)]
  if (nrow(region_dt) == 0) {
    cat("  No SNPs in region for", label, "\n")
    return(invisible(NULL))
  }
  region_dt[, log10p := -log10(P)]
  region_dt[, highlight := (P < 5e-8)]

  lead_snp_row <- region_dt[which.min(P)]

  p <- ggplot(region_dt, aes(x=BP/1e6, y=log10p)) +
    geom_point(aes(color=highlight), alpha=0.6, size=1.2) +
    geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="red") +
    geom_point(data=lead_snp_row, aes(x=BP/1e6, y=log10p),
               color="darkred", size=3, shape=18) +
    annotate("text", x=lead_snp_row$BP/1e6, y=lead_snp_row$log10p+0.5,
             label=lead_snp_row$SNP, size=3, color="darkred", fontface="bold") +
    scale_color_manual(values=c("FALSE"="steelblue","TRUE"="darkorange"), guide="none") +
    labs(
      title     = paste0("Regional Association Plot: ", label),
      subtitle  = paste0("Chr",chr,": ",round((center_bp-window_bp)/1e6,2),"-",round((center_bp+window_bp)/1e6,2)," Mb"),
      x         = paste0("Position on Chr",chr," (Mb)"),
      y         = expression(-log[10](P))
    ) +
    theme_minimal(base_size=11) +
    theme(plot.title=element_text(face="bold"))

  out_file <- file.path(FIG_DIR, paste0("regional_chr",chr,"_",label,".png"))
  ggsave(out_file, p, width=8, height=4.5, dpi=150)
  cat("  Saved:", out_file, "\n")
}

# Chr6 HLA region centred on rs9268671 (HLA-DQA1, ~32.4 Mb)
plot_region(6, 32446513, "HLA_rs9268671", window_mb=3)

# Chr6 extended HLA (top signal rs1800629, TNF locus ~31.6Mb)
plot_region(6, 31575254, "HLA_rs1800629_TNF", window_mb=2)

# Chr2 top non-HLA locus rs13019891 (STAT4, ~113 Mb)
plot_region(2, 113072292, "chr2_STAT4_rs13019891", window_mb=1)

# Chr7 IRF5 locus rs3807306
plot_region(7, 128940626, "chr7_IRF5_rs3807306", window_mb=1)

cat("\n=== Step 2 complete ===\n")
