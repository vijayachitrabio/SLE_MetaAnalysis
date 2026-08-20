library(data.table)

dir.create("data/munged", showWarnings = FALSE)

# FinnGen format: rsids, #chrom, pos, ref, alt, beta, sebeta, pval, af_alt, info
# LAVA needs: SNP, A1, A2, Z, N, P
munge_finngen <- function(file_path, trait_name, N_total) {
    cat("Harmonizing", trait_name, "...\n")
    df <- fread(file_path, select=c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"))
    
    # Filter for quality and valid RSID
    df <- df[af_alt > 0.01 & af_alt < 0.99]
    df <- df[!is.na(beta) & !is.na(sebeta) & !is.na(pval) & !is.na(rsids) & rsids != ""]
    
    # Format for LAVA: A1 must be the effect allele. FinnGen beta is for the alt allele.
    df[, Z := beta / sebeta]
    df[, N := N_total]
    setnames(df, c("rsids", "alt", "ref", "pval"), c("SNP", "A1", "A2", "P"))
    
    # Smart Deduplication
    dup_snps <- unique(df$SNP[duplicated(df$SNP)])
    if (length(dup_snps) > 0) {
        df_unique <- df[!SNP %in% dup_snps]
        df_dup <- df[SNP %in% dup_snps]
        
        # If alleles conflict, drop entirely. Else, keep lowest P.
        df_dup_resolved <- df_dup[, {
            if (uniqueN(paste(A1, A2, sep="_")) > 1) {
                .SD[0]
            } else {
                .SD[which.min(P)]
            }
        }, by = SNP]
        
        final_df <- rbind(df_unique, df_dup_resolved)[, .(SNP, A1, A2, Z, N, P)]
    } else {
        final_df <- df[, .(SNP, A1, A2, Z, N, P)]
    }
    
    out_file <- paste0("data/munged/", trait_name, "_lava_input.txt")
    fwrite(final_df, out_file, sep="\t")
    cat("Saved harmonized LAVA input for", trait_name, "to", out_file, "\n")
    return(out_file)
}

# FinnGen R12 sample sizes
ra_file <- munge_finngen("data/raw/RA_finngen_R12.gz", "RA", 315115)
ssc_file <- munge_finngen("data/raw/SSc_finngen_R12.gz", "SSc", 484260)
sjo_file <- munge_finngen("data/raw/Sjogren_finngen_R12.gz", "Sjogren", 484260)

cat("\nHarmonization complete! Ready for LAVA bivariate cross-trait analysis.\n")
