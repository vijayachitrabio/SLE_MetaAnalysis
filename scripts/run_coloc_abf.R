#!/usr/bin/env Rscript
library(coloc)
library(data.table)
library(dplyr)

cat("=========================================\n")
cat("Starting Proper Colocalization (eQTL Catalogue)\n")
cat("=========================================\n")

input_dir <- "results_extracted/coloc_inputs"
files <- list.files(input_dir, pattern="\\.tsv$", full.names=TRUE)

if (length(files) == 0) {
  stop("No input files found in ", input_dir)
}

coloc_results <- list()

for (f in files) {
  basename <- basename(f)
  parts <- strsplit(sub("\\.tsv$", "", basename), "_")[[1]]
  locus <- parts[1]
  tissue <- paste(parts[-1], collapse="_")
  
  cat(sprintf("\nProcessing %s for %s\n", locus, tissue))
  df <- fread(f)
  if (nrow(df) == 0) next
  
  # Group by gene
  genes <- unique(df$gene)
  for (g in genes) {
    g_df <- df[gene == g]
    if (nrow(g_df) < 50) next
    
    # GWAS dataset
    dataset1 <- list(
      beta = g_df$gwas_beta,
      varbeta = g_df$gwas_se^2,
      snp = g_df$snp,
      position = g_df$position,
      type = "cc",
      N = g_df$gwas_N[1],
      s = 0.0137 # 5,342 cases / 388,655 total
    )
    
    # eQTL dataset (letting coloc estimate sdY from MAF)
    dataset2 <- list(
      beta = g_df$eqtl_beta,
      varbeta = g_df$eqtl_se^2,
      snp = g_df$snp,
      position = g_df$position,
      type = "quant",
      N = g_df$eqtl_N[1],
      MAF = g_df$eqtl_maf
    )
    
    # Run COLOC
    res <- tryCatch(coloc.abf(dataset1, dataset2), error=function(e) NULL)
    
    if (!is.null(res)) {
      pp4 <- res$summary["PP.H4.abf"]
      max_p <- max(g_df$eqtl_pval, na.rm=TRUE)
      num_snps <- nrow(g_df)
      
      cat(sprintf("  -> %s: PP4 = %.3f (SNPs: %d, max p-val: %.3f)\n", g, pp4, num_snps, max_p))
      
      coloc_results[[length(coloc_results)+1]] <- data.frame(
        Locus = locus, 
        Gene = g, 
        Tissue = tissue, 
        PP4 = pp4,
        SNPs = num_snps,
        Max_eQTL_P = max_p
      )
    }
  }
}

if (length(coloc_results) > 0) {
  final_df <- bind_rows(coloc_results)
  fwrite(final_df, "results/coloc_eqtl_catalogue_summary.tsv", sep="\t")
  cat("\n=== Done! Final colocalization results saved to results/coloc_eqtl_catalogue_summary.tsv ===\n")
} else {
  cat("\nNo colocalization results generated.\n")
}
