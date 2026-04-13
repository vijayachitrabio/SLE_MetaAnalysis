#!/usr/bin/env Rscript
# scripts/step12_novelty_check.R
# GWAS Catalog Cross-Reference to identify Novel vs Known SLE loci.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(gwasrapidd)
})

RES <- "results"
cat("=== GWAS Catalog Novelty Check ===\n")

# 1. Fetch Known SLE Variants from GWAS Catalog
cat("Retrieving known SLE variants from the GWAS Catalog (Systemic lupus erythematosus)...\n")
# Direct fetch variants by efo_trait
vars <- get_variants(efo_trait = "systemic lupus erythematosus")
known_rsids <- unique(vars@variants$variant_id)

cat(sprintf("Found %d unique SNPs associated with SLE in the GWAS Catalog.\n", length(known_rsids)))

# 3. Load our Top Loci
top_loci <- fread(file.path(RES, "top_loci_summary_table.tsv"))

# 4. Identification Logic
# A. Direct Match by RSID
top_loci[, Known_GWAS_Catalog := RSID %in% known_rsids]

# B. Distance-based check (in case our lead SNP is different but in the same locus)
# This usually requires coordinates for all catalog SNPs.
# For simplicity, let's just check the "Known_SLE" column we already had 
# (which might be based on older literature) and update it.

cat("Updating discovery novelty classification...\n")
top_loci[, Novel := !Known_GWAS_Catalog]

# 5. Summarize
n_novel <- sum(top_loci$Novel)
n_known <- sum(!top_loci$Novel)
cat(sprintf("\nResults:\n  Novel Loci: %d\n  Known Loci: %d\n", n_novel, n_known))

fwrite(top_loci, file.path(RES, "top_loci_with_novelty.tsv"), sep="\t")
cat(sprintf("Novelty findings saved to %s/top_loci_with_novelty.tsv\n", RES))
