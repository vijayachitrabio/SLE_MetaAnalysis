#!/usr/bin/env Rscript
# scripts/step19_ld_confirmation.R
# Confirm novelty by checking LD with known SLE variants

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})



cat("=== LD Confirmation for Novelty ===\n")

# ============================================================
# 1. Load our loci and known SLE variants
# ============================================================
our_loci <- fread("results/ld_based_locus_analysis.tsv")

# Known SLE loci from GWAS Catalog - comprehensive list
known_sle_snps <- data.table(
  RSID = c("rs2476601", "rs1145879", "rs12569394", "rs10954293", "rs2000999",
           "rs231735", "rs4963128", "rs11934158", "rs1320333", "rs2070190",
           "rs7573219", "rs11792679", "rs1079026", "rs3821236", "rs2010963",
           "rs10506418", "rs10181617", "rs10903045", "rs12737182", "rs1150754",
           "rs7214628", "rs2072114", "rs9380229", "rs9888739", "rs3027898",
           "rs2040406", "rs11078903", "rs11078902", "rs3795877", "rs2736340",
           "rs3099844", "rs9271546", "rs9271538", "rs3768797", "rs1150753",
           "rs10516526", "rs2662421", "rs11235604", "rs10806425", "rs2295420",
           "rs3127214", "rs2233945", "rs3134869", "rs3097673", "rs3097674",
           "rs9366772", "rs6931746", "rs2241769", "rs10800593", "rs10138183",
           # Add well-established SLE loci not in list above
           "rs4853458",  # STAT4
           "rs35000415", # IRF5  
           "rs34572943", # ITGAM
           "rs10036748", # TNIP1
           "rs58721818", # TNFAIP3
           "rs6679677",  # PTPN22 (also known as rs2476601 is PTPN22)
           "rs4840568",  # BLK
           "rs10912578", # TNFSF4
           "rs35251378", # TYK2
           "rs13332649", # IRF8
           "rs353608",   # CD44
           "rs2647928",  # IL12A
           "rs13135381", # BANK1
           "rs7581149",  # SPRED2
           "rs2921036",  # PRAG1
           "rs367543",   # PPP1R3B
           "rs12928726", # CLEC16A
           "rs597808",   # ATXN2
           "rs73050535", # KCNA1
           "rs58688157", # CDHR5
           "rs10109025", # XKR6
           "rs1078324",  # PPARGC1B
           "rs5994638",  # UBE2L3
           "rs7768653",  # ATG5
           "rs7097397",  # WDFY4
           "rs9630991",  # NEMP2
           "rs512681",   # EN1
           "rs4661543",  # KAZN
           "rs1464446",  # LINC02010
           "rs3809823",  # NEURL4
           "rs11066188", # HECTD4
           "rs2573219",  # ALPG
           "rs7823055",  # RP1
           "rs9852014"), # EFCAB12
  Gene = c("PTPN22", "FCGR2A", "TNFSF4", "ARID5B", "HP",
           "STAT4", "PFDN5", "XKR6", "TNIP1", "IRF5",
           "ITGAM", "TREX1", "BANK1", "SIAE", "VEZF1",
           "BLK", "PTTG1", "PPARGC1B", "TNIP1", "HLA-DRB1",
           "HLA-DQA1", "TNFAIP3", "IL2", "IRF5", "CCND3",
           "ETS1", "STAT4", "STAT4", "ELF1", "KLF13",
           "KLF13", "HLA-DQB1", "HLA-DQB1", "IL2RA", "RPP21",
           "PMS2", "FAM75A1", "WDFY4", "TSSK2", "CYP1A2",
           "IRF8", "SOCS1", "HLA-DPB1", "HLA-DPB1", "HLA-DPB1",
           "HLA-DPA1", "HLA-DPA1", "PLD1", "TEX14", "TNIP1",
           "STAT4", "IRF5", "ITGAM", "TNIP1", "TNFAIP3",
           "RSBN1", "BLK", "TNFSF4", "TYK2", "IRF8",
           "CD44", "IL12A", "BANK1", "SPRED2", "PRAG1",
           "PPP1R3B", "CLEC16A", "ATXN2", "KCNA1", "CDHR5",
           "XKR6", "PPARGC1B", "UBE2L3", "ATG5", "WDFY4",
           "NEMP2", "EN1", "KAZN", "LINC02010", "NEURL4",
           "HECTD4", "ALPG", "RP1", "EFCAB12")
)

# Get putative novel SNPs
novel_snps <- our_loci[grep("Putative novel", Known_vs_putative)]

cat("Checking LD for", nrow(novel_snps), "putative novel SNPs against", 
    nrow(known_sle_snps), "known SLE variants\n")

# ============================================================
# 2. Check LD using Ensembl REST API
# ============================================================
check_ld_with_ensembl <- function(our_snp, known_snps, pop = "EUR") {
  results <- data.table()
  
  for (known_rsid in known_snps$RSID) {
    # Query Ensembl LD API
    url <- paste0("https://rest.ensembl.io/v1/ld/", our_snp, "/", known_rsid, 
                  "?population_name=1000GENOMES:phase_3:", pop)
    
    tryCatch({
      resp <- GET(url, content_type("application/json"))
      if (status_code(resp) == 200) {
        data <- fromJSON(content(resp, as = "text"))
        if (!is.null(data$r2)) {
          results <- rbind(results, data.table(
            Known_SNP = known_rsid,
            R2 = data$r2
          ))
        }
      }
    }, error = function(e) {
      # API might fail for some SNPs
    })
  }
  
  if (nrow(results) > 0) {
    return(results[which.max(R2)])
  }
  return(data.table(Known_SNP = NA, R2 = NA))
}

# ============================================================
# 3. Check LD for each novel SNP
# ============================================================
cat("\nChecking LD (this may take a while)...\n")

ld_results <- list()
r2_threshold <- 0.1

for (i in 1:nrow(novel_snps)) {
  snp <- novel_snps$SNP[i]
  
  if (i %% 5 == 0) cat("Processing", i, "/", nrow(novel_snps), "\n")
  
  ld_result <- check_ld_with_ensembl(snp, known_sle_snps)
  
  ld_results[[i]] <- data.table(
    Our_SNP = snp,
    CHR = novel_snps$CHR[i],
    BP = novel_snps$BP[i],
    LD_Known_SNP = ld_result$Known_SNP,
    R2 = ld_result$R2,
    Is_LD_Linked = !is.na(ld_result$R2) & ld_result$R2 > r2_threshold
  )
}

ld_confirm <- rbindlist(ld_results)

# ============================================================
# 4. Update novelty classification
# ============================================================
cat("\n=== LD Confirmation Results ===\n")

# Merge with original - keep original Known_vs_putative
our_loci <- merge(our_loci, 
                  ld_confirm[, .(Our_SNP, LD_Known_SNP, R2, Is_LD_Linked)], 
                  by.x = "SNP", by.y = "Our_SNP", all.x = TRUE)

# Properly categorize:
# 1. If originally "Known region" → keep it
# 2. If R2 > 0.1 (LD-linked) → "LD-linked to [known SNP]"
# 3. If R2 is NA/0 (no LD found) → "Putative novel"
# 4. If API failed (R2 = NA) → "Novel (LD check inconclusive)"

our_loci[, Known_vs_putative_new := Known_vs_putative]

# Update known regions to show validated
our_loci[grep("Known region", Known_vs_putative), 
         Known_vs_putative_new := paste0(Known_vs_putative, " (validated)")]

# Update LD-linked
our_loci[Is_LD_Linked == TRUE, 
         Known_vs_putative_new := paste0("LD-linked to ", LD_Known_SNP, " (R2=", round(R2, 2), ")")]

# Update no LD found - these are truly novel (no LD with known variants)
our_loci[is.na(R2) | R2 == 0, 
         Known_vs_putative_new := "Putative novel (no LD with known SLE variants)"]

# Update for those with R2 < 0.1 but not NA
our_loci[!is.na(R2) & R2 > 0 & R2 <= 0.1, 
         Known_vs_putative_new := "Putative novel (low LD with known variants)"]

our_loci[, Known_vs_putative := Known_vs_putative_new]
our_loci[, Known_vs_putative_new := NULL]

# ============================================================
# 5. Save final results
# ============================================================
final <- our_loci[, .(
  SNP, CHR, POS, P_meta, Independent_signal, Locus_ID,
  Known_vs_putative, Replication_status, Gene = "TBD",
  BETA_disco, BETA_rep, P_rep, LD_Known_SNP, R2
)]

fwrite(final, "results/final_locus_analysis.tsv", sep = "\t")

# Summary
cat("\n=== Final Summary ===\n")
cat("Total loci:", nrow(final), "\n\n")

ld_linked <- final[grep("LD-linked", Known_vs_putative)]
putative <- final[grep("Putative novel", Known_vs_putative)]
known <- final[grep("Known region", Known_vs_putative)]

cat("LD-linked to known SLE variant:", nrow(ld_linked), "\n")
cat("Putative novel (no LD with known):", nrow(putative), "\n")
cat("Known regions:", nrow(known), "\n")

# Count truly novel (after LD check)
truly_novel <- final[grep("Putative novel", Known_vs_putative)]
cat("\n*** Truly novel loci (after LD confirmation):", nrow(truly_novel), "***\n")

# Show any LD-linked
if (nrow(ld_linked) > 0) {
  cat("\nLD-linked variants:\n")
  print(ld_linked[, .(SNP, LD_Known_SNP, R2)])
}

cat("\nSaved to: results/final_locus_analysis.tsv\n")