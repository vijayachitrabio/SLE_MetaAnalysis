#!/usr/bin/env Rscript
# scripts/step24_coloc_analysis.R
# Functional Audit: Colocalization (COLOC) for 14 High-Confidence Loci

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})

cat("=========================================\n")
cat("Starting Colocalization Analysis (COLOC)\n")
cat("=========================================\n")

# 1. Load High-Confidence Loci
lava_loci <- fread("results/final_lava_consolidated_loci.tsv")[Final_Assessment == "HIGH CONFIDENCE"]
gwas_full <- fread("results/discovery_meta_results.tsv")
# Reference sample size for GWAS (Meta): Corrected to audited Discovery counts
N_GWAS <- 388655 
S_GWAS <- 0.0137 # 5,342 cases / 388,655 total

# 2. API Helper Functions
get_ensembl_genes <- function(chr, start, end) {
  # Query Ensembl REST API for genes in range
  url <- sprintf("https://rest.ensembl.org/overlap/region/human/%s:%d-%d?feature=gene;content-type=application/json", chr, start, end)
  res <- tryCatch(GET(url), error=function(e) NULL)
  if (!is.null(res) && status_code(res) == 200) {
    genes <- fromJSON(content(res, as = "text"))
    if (length(genes) > 0) return(unique(genes$external_name))
  }
  return(NULL)
}

get_gtex_eqtl <- function(gene_symbol, tissue="Whole_Blood") {
  # Get Gencode ID first
  url_gene <- "https://gtexportal.org/api/v2/reference/gene"
  res_gene <- tryCatch(GET(url_gene, query = list(geneId = gene_symbol, datasetId = "gtex_v10")), error=function(e) NULL)
  if (is.null(res_gene) || status_code(res_gene) != 200) return(NULL)
  gene_data <- fromJSON(content(res_gene, as = "text"))$data
  if (length(gene_data) == 0) return(NULL)
  gencode_id <- gene_data$gencodeId[1]
  
  # Fetch Associations
  url_assoc <- "https://gtexportal.org/api/v2/association/singleTissueEqtl"
  res_assoc <- tryCatch(GET(url_assoc, query = list(
    gencodeId = gencode_id,
    tissueSiteDetailId = tissue,
    datasetId = "gtex_v10"
  )), error=function(e) NULL)
  
  if (!is.null(res_assoc) && status_code(res_assoc) == 200) {
    return(fromJSON(content(res_assoc, as = "text"))$data)
  }
  return(NULL)
}

# 3. Main Analysis Loop
coloc_all_results <- list()

for (i in 1:nrow(lava_loci)) {
  rsid <- lava_loci$RSID[i]
  chr <- lava_loci$CHR[i]
  pos <- lava_loci$BP[i]
  
  cat(sprintf("\nProcessing Locus: %s (Chr%d:%d)\n", rsid, chr, pos))
  
  # A. Define Window
  start <- pos - 250000
  end <- pos + 250000
  
  # B. Get Genes in Window
  genes <- get_ensembl_genes(chr, start, end)
  if (is.null(genes) || length(genes) == 0) {
    # Fallback to the one mapped gene if available
    genes <- lava_loci$Gene[i]
  }
  
  # Ensure genes is not NULL and remove NA/empty
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) {
    cat("  No genes identified in window.\n")
    next
  }
  
  # C. GWAS Window data
  gwas_subset <- gwas_full[CHR == chr & BP >= start & BP <= end & !is.na(P_meta) & !is.na(RSID)]
  if (nrow(gwas_subset) < 50) next
  
  d1 <- list(
    beta = gwas_subset$BETA_meta,
    varbeta = gwas_subset$SE_meta^2,
    snp = gwas_subset$RSID,
    position = gwas_subset$BP,
    type = "cc",
    N = N_GWAS,
    s = S_GWAS
  )

  # D. Loop through Genes and Tissues
  for (gene in unique(genes)) {
    # Try Whole Blood and Spleen
    for (tis in c("Whole_Blood", "Spleen")) {
      cat(sprintf("  Querying %s in %s...\n", gene, tis))
      eqtl <- tryCatch(get_gtex_eqtl(gene, tis), error=function(e) NULL)
      
      if (is.null(eqtl) || !is.data.frame(eqtl)) {
        cat("    No valid eQTL data frame returned.\n")
        next
      }
      if (nrow(eqtl) < 50) {
        cat(sprintf("    Too few associations (%d) for COLOC.\n", nrow(eqtl)))
        next
      }
      
      # Intersect
      common <- intersect(d1$snp, eqtl$snpId)
      if (length(common) < 50) next
      
      idx1 <- match(common, d1$snp)
      idx2 <- match(common, eqtl$snpId)
      
      dataset1 <- list(
        beta = d1$beta[idx1],
        varbeta = d1$varbeta[idx1],
        snp = common,
        position = d1$position[idx1],
        type = "cc", N = N_GWAS, s = d1$s
      )
      # Approximate SE from NES and P-value
      pvals <- as.numeric(eqtl$pValue[idx2])
      pvals[pvals < 1e-300] <- 1e-300 # Cap small P-values
      pvals[pvals > 0.999] <- 0.999   # Cap large P-values
      
      dataset2 <- list(
        beta = as.numeric(eqtl$nes[idx2]),
        varbeta = (as.numeric(eqtl$nes[idx2]) / abs(qnorm(pvals/2)))^2,
        snp = common,
        position = as.numeric(eqtl$pos[idx2]),
        type = "quant",
        N = 670, # GTEx v8
        sdY = 1 # Standard for normalized GTEx data
      )
      
      # Run COLOC
      res <- coloc.abf(dataset1, dataset2)
      pp4 <- res$summary["PP.H4.abf"]
      cat(sprintf("    -> PP4 = %.3f\n", pp4))
      
      if (pp4 > 0.8) {
        coloc_all_results[[length(coloc_all_results)+1]] <- data.frame(
          Locus = rsid, Gene = gene, Tissue = tis, PP4 = pp4
        )
      }
    }
  }
  Sys.sleep(1) # API kindness
}

# 4. Save
if (length(coloc_all_results) > 0) {
  final_df <- bind_rows(coloc_all_results)
  fwrite(final_df, "results/coloc_results_summary.tsv", sep="\t")
  cat("\n=== Done! Final colocalization results saved to results/coloc_results_summary.tsv ===\n")
} else {
  cat("\nNo significant colocalization (PP4 > 0.5) found.\n")
}
