#!/usr/bin/env Rscript
# scripts/step28_immune_epigenetic_overlap.R
# Exploratory Immune Chromatin-State Overlap (Roadmap Epigenomics)

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(rtracklayer)
})

cat("=========================================\n")
cat("Starting Exploratory Epigenetic Overlap\n")
cat("=========================================\n")

# 1. Download chain file for hg38 -> hg19 liftover if it doesn't exist
chain_file <- "reference_data/hg38ToHg19.over.chain"
if (!file.exists(chain_file)) {
    cat("Downloading UCSC hg38ToHg19 liftOver chain...\n")
    chain_url <- "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz"
    download.file(chain_url, paste0(chain_file, ".gz"), quiet = TRUE)
    system(paste("gunzip", paste0(chain_file, ".gz")))
}
chain <- import.chain(chain_file)

# 2. Load 47 top loci from master results
cat("Loading 47 SLE discovery loci...\n")
master_res <- fread("results_extracted/master_results_table.tsv")

# Identify MHC loci (chr6: 25Mb-35Mb)
master_res[, Is_MHC := (CHR == 6 & BP >= 25000000 & BP <= 35000000)]
cat(sprintf("Found %d non-MHC loci and %d MHC loci.\n", sum(!master_res$Is_MHC), sum(master_res$Is_MHC)))

# Create GRanges (hg38)
gr_hg38 <- GRanges(
    seqnames = paste0("chr", master_res$CHR),
    ranges = IRanges(start = master_res$BP, end = master_res$BP),
    LOC = master_res$LAVA_LOC,
    RSID = master_res$RSID,
    Is_MHC = master_res$Is_MHC
)

# 3. LiftOver to hg19
cat("Lifting coordinates to hg19...\n")
gr_hg19_list <- liftOver(gr_hg38, chain)
# For SNPs, keep those that mapped to exactly 1 position
mapped <- lengths(gr_hg19_list) == 1
cat(sprintf("Successfully lifted %d out of %d loci to hg19.\n", sum(mapped), length(mapped)))

gr_hg19_exact <- unlist(gr_hg19_list[mapped])

# Create ±10kb proximal regulatory window
gr_hg19_window <- GRanges(
    seqnames = seqnames(gr_hg19_exact),
    ranges = IRanges(start = start(gr_hg19_exact) - 10000, end = end(gr_hg19_exact) + 10000),
    LOC = gr_hg19_exact$LOC,
    RSID = gr_hg19_exact$RSID,
    Is_MHC = gr_hg19_exact$Is_MHC
)

# 4. Roadmap Data (E062=PBMC, E032=B cells, E034=T cells, E029=Monocytes)
cell_types <- c("E062", "E032", "E034", "E029")
cell_names <- c("E062_PBMC", "E032_B_cell", "E034_T_cell", "E029_Monocyte")
roadmap_url_base <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"

results_list <- list()

# Active states: 1_TssA, 2_TssAFlnk, 6_EnhG, 7_Enh 
target_states <- c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh")

cat("Downloading and overlapping with Roadmap 15-state Epigenomes...\n")

overlap_results <- data.table(
    LOC = gr_hg19_exact$LOC,
    RSID = gr_hg19_exact$RSID,
    Is_MHC = gr_hg19_exact$Is_MHC
)

for (i in 1:length(cell_types)) {
    eid <- cell_types[i]
    cname <- cell_names[i]
    bed_file <- paste0("reference_data/", eid, "_15_coreMarks_mnemonics.bed")
    
    if (!file.exists(bed_file)) {
        url <- paste0(roadmap_url_base, eid, "_15_coreMarks_mnemonics.bed.gz")
        # Use curl to bypass SSL issues natively
        system(sprintf("curl -kLs %s -o %s.gz", url, bed_file))
        system(paste("gunzip", paste0(bed_file, ".gz")))
    }
    
    # Read BED (chr, start, end, state)
    bed <- fread(bed_file, col.names = c("chr", "start", "end", "state"))
    # Filter for active states
    bed_active <- bed[state %in% target_states]
    
    gr_bed <- GRanges(seqnames = bed_active$chr, ranges = IRanges(start = bed_active$start, end = bed_active$end))
    
    # Overlap Exact
    hits_exact <- findOverlaps(gr_hg19_exact, gr_bed)
    loc_has_exact <- rep(FALSE, length(gr_hg19_exact))
    loc_has_exact[queryHits(hits_exact)] <- TRUE
    
    # Overlap Window (±10kb)
    hits_window <- findOverlaps(gr_hg19_window, gr_bed)
    loc_has_window <- rep(FALSE, length(gr_hg19_window))
    loc_has_window[queryHits(hits_window)] <- TRUE
    
    overlap_results[[paste0(cname, "_Exact")]] <- loc_has_exact
    overlap_results[[paste0(cname, "_Window10kb")]] <- loc_has_window
    
    cat(sprintf("  [%s] Exact overlap: %d | Window overlap: %d\n", cname, sum(loc_has_exact), sum(loc_has_window)))
}

# Aggregate: Does it overlap ANY of the 4 immune cells?
overlap_results[, Any_Immune_Exact := (E062_PBMC_Exact | E032_B_cell_Exact | E034_T_cell_Exact | E029_Monocyte_Exact)]
overlap_results[, Any_Immune_Window10kb := (E062_PBMC_Window10kb | E032_B_cell_Window10kb | E034_T_cell_Window10kb | E029_Monocyte_Window10kb)]

fwrite(overlap_results, "results_extracted/exploratory_epigenetic_overlap.tsv", sep="\t")

# Output Summary for Non-MHC
non_mhc <- overlap_results[Is_MHC == FALSE]
total_non_mhc <- nrow(non_mhc)

cat("\n=========================================\n")
cat("SUMMARY FOR NON-MHC LOCI (n =", total_non_mhc, ")\n")
cat("=========================================\n")
cat(sprintf("Exact Lead-SNP Overlap (Any Immune Cell): %d (%.1f%%)\n", sum(non_mhc$Any_Immune_Exact), 100 * mean(non_mhc$Any_Immune_Exact)))
cat(sprintf("Proximal Window ±10kb Overlap (Any Immune Cell): %d (%.1f%%)\n", sum(non_mhc$Any_Immune_Window10kb), 100 * mean(non_mhc$Any_Immune_Window10kb)))

cat("\nDone! Results saved to results_extracted/exploratory_epigenetic_overlap.tsv\n")
