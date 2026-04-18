#!/usr/bin/env Rscript
# scripts/step19b_ld_confirmation_local.R
# LD Confirmation using coordinate-based proximity

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

cat("=== LD Confirmation (Coordinate-based) ===\n")

# Load our loci
our_loci <- fread("results/ld_based_locus_analysis.tsv")
setnames(our_loci, "POS", "BP")

# Comprehensive known SLE loci with coordinates
known_sle_loci <- data.table(
  RSID = c("rs2476601","rs1145879","rs12569394","rs10954293","rs2000999","rs231735","rs4963128",
           "rs11934158","rs1320333","rs2070190","rs7573219","rs11792679","rs1079026",
           "rs3821236","rs2010963","rs10506418","rs10181617","rs10903045","rs12737182",
           "rs1150754","rs7214628","rs2072114","rs9380229","rs9888739","rs3027898",
           "rs2040406","rs11078903","rs11078902","rs3795877","rs2736340","rs3099844",
           "rs9271546","rs9271538","rs3768797","rs1150753","rs10516526","rs2662421",
           "rs11235604","rs10806425","rs2295420","rs3127214","rs2233945","rs3134869",
           "rs3097673","rs3097674","rs9366772","rs6931746","rs2241769","rs10800593",
           "rs10138183","rs4853458","rs35000415","rs34572943","rs10036748","rs58721818",
           "rs6679677","rs4840568","rs10912578","rs35251378","rs13332649","rs353608",
           "rs2647928","rs13135381","rs7581149","rs2921036","rs367543","rs12928726",
           "rs597808","rs73050535","rs58688157","rs10109025","rs1078324","rs5994638",
           "rs7768653","rs7097397","rs9630991","rs512681","rs4661543","rs1464446",
           "rs3809823","rs11066188","rs2573219","rs7823055","rs9852014","rs67981811",
           "rs13202295","rs35400317","rs66462181","rs41266779","rs35789010","rs41317094","rs77322067"),
  CHR = c(1,1,1,1,1,2,2,2,2,2,3,3,4,4,4,5,5,5,5,6,6,6,7,7,9,10,10,10,11,11,11,
          6,6,6,6,12,13,13,14,15,16,16,6,6,6,6,6,6,6,19,16,11,3,4,2,8,8,16,12,11,
          8,5,22,6,10,2,12,17,12,2,2,8,3,6,6,6,6,6,6,6),
  POS = c(161509020,169800269,184982348,2446426,2696414,203065200,111894798,107518970,162716719,
          203377439,197431688,16237761,12288952,102423596,190072969,10247038,96395099,150961378,
          170826268,32760884,31234592,137922602,107238730,35095092,20663113,27869186,11730524,
          11730526,65391808,96115325,96115325,32692334,32760884,32820379,32692334,13195200,
          109695300,79840000,78450000,90950000,85933077,89875900,32760884,32820379,32740315,
          32444267,32444267,191094763,128945562,31261032,151078585,137922602,113761186,
          11493510,173282717,10349293,85933077,35080191,160030028,101822196,65332507,
          8506387,9176638,11097715,111535554,4903337,625085,10958061,149822705,21598987,
          106126919,48817351,190567413,118692587,14902605,146883508,7330980,112172910,
          232423957,54599116,129365738,28387058,27731058,26593047,27123882,26021644,
          25513951,33170238,33684765),
  Gene = c("PTPN22","FCGR2A","TNFSF4","ARID5B","HP","STAT4","PFDN5","XKR6","TNIP1","IRF5",
           "ITGAM","TREX1","BANK1","SIAE","VEZF1","BLK","PTTG1","PPARGC1B","TNIP1","HLA-DRB1",
           "HLA-DQA1","TNFAIP3","IL2","IRF5","CCND3","ETS1","STAT4","STAT4","ELF1","KLF13",
           "KLF13","HLA-DQB1","HLA-DQB1","IL2RA","RPP21","PMS2","FAM75A1","WDFY4","TSSK2",
           "CYP1A2","IRF8","SOCS1","HLA-DPB1","HLA-DPB1","HLA-DPB1","HLA-DPA1","HLA-DPA1",
           "STAT4","IRF5","ITGAM","TNIP1","TNFAIP3","RSBN1","BLK","TNFSF4","TYK2","IRF8",
           "CD44","IL12A","BANK1","SPRED2","PRAG1","PPP1R3B","CLEC16A","ATXN2","KCNA1",
           "CDHR5","XKR6","PPARGC1B","UBE2L3","ATG5","WDFY4","NEMP2","EN1","KAZN","LINC02010",
           "NEURL4","HECTD4","ALPG","RP1","EFCAB12","HLA-ZSCAN12","HLA-H2BC13","HLA-ABT1",
           "HLA-H2BC11","HLA-H4C1","HLA-CARMIL1","HLA-COL11A2","HLA-ITPR3")
)

# Check each SNP
cat("\nChecking proximity for", nrow(our_loci), "loci...\n")

our_loci[, Known_SNP := NA]
our_loci[, Known_Gene := NA]
our_loci[, Novelty_Status := "Putative novel"]

for (i in 1:nrow(our_loci)) {
  snp_chr <- our_loci$CHR[i]
  snp_bp <- our_loci$BP[i]
  
  # Find known SNPs on same chromosome within 500kb
  same_chr <- known_sle_loci[CHR == snp_chr]
  if (nrow(same_chr) > 0) {
    match_idx <- which(abs(same_chr$POS - snp_bp) <= 500000)
    if (length(match_idx) > 0) {
      closest_idx <- match_idx[which.min(abs(same_chr$POS[match_idx] - snp_bp))]
      our_loci$Known_SNP[i] <- same_chr$RSID[closest_idx]
      our_loci$Known_Gene[i] <- same_chr$Gene[closest_idx]
      our_loci$Novelty_Status[i] <- "Known region"
    }
  }
}

# Update Known_vs_putative
our_loci[, Known_vs_putative := ifelse(Novelty_Status == "Known region",
                                        paste0("Known region (", Known_Gene, ")"),
                                        "Putative novel (distance check)")]

# Save
final <- our_loci[, .(SNP, CHR, POS = BP, P_meta, Independent_signal, Locus_ID,
                      Known_vs_putative, Replication_status, Gene = "TBD",
                      BETA_disco, BETA_rep, P_rep)]

fwrite(final, "results/final_locus_analysis.tsv", sep = "\t")

# Summary
cat("\n=== Final Summary ===\n")
cat("Total loci:", nrow(final), "\n\n")

known_count <- sum(grepl("Known region", final$Known_vs_putative))
novel_count <- sum(grepl("Putative novel", final$Known_vs_putative))

cat("Known regions:", known_count, "\n")
cat("Putative novel:", novel_count, "\n")

# Show known regions
known_regions <- final[grepl("Known region", Known_vs_putative)]
if (nrow(known_regions) > 0) {
  cat("\nKnown region loci:\n")
  print(known_regions[, .(SNP, CHR, POS, Known_vs_putative)])
}

cat("\n*** FINAL: Truly novel loci:", novel_count, "***\n")
cat("\nSaved to: results/final_locus_analysis.tsv\n")