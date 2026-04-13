library(data.table)
library(dplyr)
library(ggplot2)
library(gprofiler2)

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

message("Loading top annotated loci...")
top_loci <- fread("results/top_loci_summary_table.tsv")

# Extract unique genes from the Gene column
# Some are like "HLA / ZSCAN12", extract the actual gene
genes <- top_loci$Gene
genes <- sapply(strsplit(genes, " / "), function(x) x[length(x)])
genes <- unique(genes[genes != "TBD" & genes != "Intergenic"])

message(paste("Found", length(genes), "unique mapped genes for enrichment."))

message("Running gProfiler enrichment...")
gostres <- gost(query = genes, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO:BP", "REAC", "KEGG"))

if (is.null(gostres)) {
  message("No significant pathways found.")
} else {
  message("Plotting pathway enrichment results...")
  p_gost <- gostplot(gostres, capped = TRUE, interactive = FALSE)
  ggsave("figures/pathway_enrichment.png", p_gost, width = 10, height = 6, dpi=300)
  
  message("Saving enrichment table...")
  enrich_df <- gostres$result %>% select(source, term_id, term_name, p_value, intersection_size) %>% arrange(p_value)
  fwrite(enrich_df, "results/pathway_enrichment_results.tsv", sep="\t")
}

message("Done with step 6.")
