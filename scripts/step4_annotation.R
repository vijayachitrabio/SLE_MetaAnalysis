library(vroom)
library(dplyr)
library(biomaRt)

# Set working directory to project root
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

message("Loading top loci summary table...")
path_table <- "results/top_loci_summary_table.tsv"
final_table <- vroom(path_table)

message("Connecting to Ensembl biomaRt (hg38)...")
# Use grch38 ensembl with mirror to avoid connection issues
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Function to get nearest gene within +/- 250kb
get_nearest_gene <- function(chr, bp) {
  # Add slightly larger window for boundary genes
  window <- 250000 
  start_pos <- max(1, bp - window)
  end_pos <- bp + window
  
  query <- getBM(
    attributes = c('external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'gene_biotype'),
    filters = c('chromosome_name', 'start', 'end'),
    values = list(chr, start_pos, end_pos),
    mart = ensembl
  )
  
  # Filter for protein coding mostly, or at least have a gene name
  query <- query %>% filter(external_gene_name != "")
  
  if (nrow(query) == 0) {
    return("Intergenic")
  }
  
  # Calculate distance
  query <- query %>%
    mutate(
      dist_start = abs(start_position - bp),
      dist_end = abs(end_position - bp),
      min_dist = pmin(dist_start, dist_end)
    ) %>%
    # If the variant is inside the gene, dist is 0
    mutate(min_dist = ifelse(bp >= start_position & bp <= end_position, 0, min_dist)) %>%
    arrange(min_dist)
  
  # Select top nearest 1-2 genes (protein coding preferred)
  pc_query <- query %>% filter(gene_biotype == "protein_coding")
  if(nrow(pc_query) > 0) {
    # If there's a protein coding gene, return it
    top_gene <- pc_query$external_gene_name[1]
  } else {
    top_gene <- query$external_gene_name[1]
  }
  
  return(top_gene)
}

message("Annotating genes for 57 loci...")
# We only want to annotate TBD or empty genes, but doing it for all ensures consistency.
# Special handling for HLA which spans many genes
for (i in 1:nrow(final_table)) {
  if (final_table$CHR[i] == 6 && final_table$BP[i] > 25000000 && final_table$BP[i] < 35000000) {
    final_table$Gene[i] <- paste0("HLA / ", get_nearest_gene(final_table$CHR[i], final_table$BP[i]))
  } else {
    gene <- get_nearest_gene(final_table$CHR[i], final_table$BP[i])
    final_table$Gene[i] <- gene
  }
}

# Add region info simply by cytoband format loosely (e.g. 1q21) - we can skip if not strictly needed or infer.
# We'll just leave region as genomic coord or simple chromosome
final_table$Region <- paste0(final_table$CHR, "p/q")

message("Saving annotated table...")
vroom_write(final_table, path_table)
message("Done with annotation.")
