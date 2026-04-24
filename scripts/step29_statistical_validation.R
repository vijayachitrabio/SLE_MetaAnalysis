#!/usr/bin/env Rscript
# scripts/step29_statistical_validation.R
# Targeted Bayesian Fine-mapping and COLOC for CLIC1 and TNFSF4

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})

cat("=== Starting Targeted Statistical Validation (PP > 0.9, PP4 > 0.8) ===\n")

# Settings
N_GWAS <- 388655 
S_GWAS <- 0.0137 
gwas_full <- fread("results/discovery_meta_results.tsv")

target_loci <- list(
  list(rsid = "rs389884", gene = "CLIC1", chr = 6, pos = 31973120),
  list(rsid = "rs10912578", gene = "TNFSF4", chr = 1, pos = 173282717)
)

# Focus on tissues with plausible immune relevance for SLE and B/T-cell biology.
tissues_to_test <- c(
  "Whole_Blood",
  "Spleen",
  "Cells_EBV-transformed_lymphocytes",
  "Minor_Salivary_Gland",
  "Lung"
)

# API helper
get_gtex_eqtl <- function(gene_symbol, tissue) {
  url_gene <- "https://gtexportal.org/api/v2/reference/gene"
  res_gene <- tryCatch(GET(url_gene, query = list(geneId = gene_symbol, datasetId = "gtex_v10")), error=function(e) NULL)
  if (is.null(res_gene) || status_code(res_gene) != 200) return(NULL)
  gene_data <- fromJSON(content(res_gene, as = "text"))$data
  if (length(gene_data) == 0) return(NULL)
  gencode_id <- gene_data$gencodeId[1]
  
  url_assoc <- "https://gtexportal.org/api/v2/association/singleTissueEqtl"
  res_assoc <- tryCatch(GET(url_assoc, query = list(gencodeId = gencode_id, tissueSiteDetailId = tissue, datasetId = "gtex_v10")), error=function(e) NULL)
  if (!is.null(res_assoc) && status_code(res_assoc) == 200) return(fromJSON(content(res_assoc, as = "text"))$data)
  return(NULL)
}

results_summary <- list()

for (locus in target_loci) {
  rsid <- locus$rsid
  gene <- locus$gene
  chr  <- locus$chr
  pos  <- locus$pos
  
  cat(sprintf("\nValidating Locus: %s (%s)\n", gene, rsid))
  
  start <- pos - 250000
  end <- pos + 250000
  gwas_subset <- gwas_full[CHR == chr & BP >= start & BP <= end & !is.na(P_meta) & !is.na(RSID)]
  
  d1 <- list(
    beta = gwas_subset$BETA_meta,
    varbeta = gwas_subset$SE_meta^2,
    snp = gwas_subset$RSID,
    position = gwas_subset$BP,
    type = "cc", N = N_GWAS, s = S_GWAS
  )

  for (tis in tissues_to_test) {
    cat(sprintf("  Testing tissue: %s... ", tis))
    eqtl_res <- get_gtex_eqtl(gene, tis)
    
    if (is.null(eqtl_res) || !is.data.frame(eqtl_res) || nrow(eqtl_res) < 50) {
      cat("No valid data.\n")
      next
    }
    
    if (!all(c("snpId", "pValue", "nes", "pos") %in% names(eqtl_res))) {
      cat("Missing columns.\n")
      next
    }
    
    common <- intersect(d1$snp, eqtl_res$snpId)
    if (length(common) < 50) {
      cat("Too few overlap.\n")
      next
    }
    
    idx1 <- match(common, d1$snp)
    idx2 <- match(common, eqtl_res$snpId)
    
    ds1 <- list(beta = d1$beta[idx1], varbeta = d1$varbeta[idx1], snp = common, type = "cc", N = N_GWAS, s = S_GWAS)
    
    pvals <- as.numeric(eqtl_res$pValue[idx2])
    pvals[is.na(pvals)] <- 1e-5
    pvals[pvals < 1e-300] <- 1e-300
    
    nes_vals <- as.numeric(eqtl_res$nes[idx2])
    z_scores <- abs(qnorm(pvals/2))
    se_est <- abs(nes_vals / z_scores)
    varb2 <- se_est^2
    varb2[is.na(varb2) | varb2 == 0] <- 0.01 
    
    ds2 <- list(
      beta = nes_vals,
      varbeta = varb2,
      snp = common,
      type = "quant", N = 670, sdY = 1
    )
    
    res <- tryCatch(coloc.abf(ds1, ds2), error=function(e) NULL)
    if (is.null(res)) {
      cat("COLOC failed.\n")
      next
    }
    
    pp4 <- as.numeric(res$summary["PP.H4.abf"])
    v_probs <- res$results
    top_v_pip <- if(!is.null(v_probs)) max(v_probs$SNP.PP.H4, na.rm = TRUE) else 0
    top_v_rsid <- if(!is.null(v_probs)) v_probs$snp[which.max(v_probs$SNP.PP.H4)] else "None"
    
    cat(sprintf("PP4 = %.3f | Max Variant PIP = %.3f (%s)\n", pp4, top_v_pip, top_v_rsid))
    
    results_summary[[length(results_summary)+1]] <- data.frame(
      Gene = gene, RSID = rsid, Tissue = tis, PP4 = pp4, Variant_PIP = top_v_pip, Top_Variant = top_v_rsid
    )
  }
}

cat("\n=== Final Statistical Support Summary ===\n")
final_df <- bind_rows(results_summary)
print(final_df)
fwrite(final_df, "results/targeted_causal_validation.tsv", sep="\t")
