#!/usr/bin/env Rscript
# scripts/step7_ppi_network.R
# Generate a high-quality Protein-Protein Interaction (PPI) network using the 57 lead loci
# and enriched pathway drivers from the new discovery analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(dplyr)
  library(httr)
  library(jsonlite)
  library(stringr)
})

RESULTS_DIR <- "results"
FIG_DIR     <- "figures"

FILE_ANNOT    <- file.path(RESULTS_DIR, "top_loci_summary_table.tsv")
FILE_PATHWAY  <- file.path(RESULTS_DIR, "pathway_enrichment_results.tsv")
OUTPUT_PNG    <- file.path(FIG_DIR, "Fig_3_PPI_network_new.png")

cat("=== Building Enhanced PPI Network ===\n")

if (!file.exists(FILE_ANNOT)) stop("Gene annotation file not found.")
if (!file.exists(FILE_PATHWAY)) stop("Pathway results file not found.")

# A. Lead SNP Genes
annot_dt <- fread(FILE_ANNOT)
annot_dt <- annot_dt[Gene != "TBD" & Gene != "Intergenic" & !is.na(Gene)]

seed_genes_df <- annot_dt %>%
  mutate(Gene = str_split(Gene, " / ") %>% sapply(function(x) x[length(x)])) %>% 
  group_by(Gene) %>%
  summarize(log10P = max(-log10(P_disco), na.rm = TRUE), .groups = "drop")

# B. Pathway Driver Genes (Extract from top significant pathways)
path_dt <- fread(FILE_PATHWAY)
driver_genes <- c()
# If gprofiler didn't output intersection string natively we will skip specific drivers
# and just use the lead genes. Let's check if 'intersection' column exists.
if("intersection" %in% colnames(path_dt)) {
  driver_genes <- path_dt %>%
    arrange(p_value) %>%
    slice_head(n = 5) %>%
    pull(intersection) %>%
    paste(collapse = ",") %>%
    str_split(",") %>%
    unlist() %>%
    unique()
}

# C. Merge and Label
all_genes_meta <- data.frame(Gene = unique(c(seed_genes_df$Gene, driver_genes))) %>%
  left_join(seed_genes_df, by = "Gene") %>%
  mutate(
    type = ifelse(!is.na(log10P), "Lead SNP Gene", "Pathway Driver"),
    log10P = ifelse(is.na(log10P), 2, log10P) # Default size for pathway genes
  )

cat("Genes identified:", nrow(all_genes_meta), "(Lead:", sum(all_genes_meta$type == "Lead SNP Gene"), ")\n")

# --- 2. Fetch Interactions from STRING API ---
cat("Fetching PPI data from STRING API (Species 9606, human)...\n")
string_api_url <- "https://string-db.org/api/json/network"
res <- GET(string_api_url, query = list(
  identifiers = paste(all_genes_meta$Gene, collapse = "\r"),
  species = 9606,
  required_score = 400, # Medium confidence to allow more connections for 57 genes
  caller_identity = "sle_meta_analysis_v2"
))

if (status_code(res) != 200) {
  stop("Failed to fetch data from STRING API. Status: ", status_code(res))
}

interactions <- fromJSON(content(res, as = "text"))
if (length(interactions) == 0) {
  stop("No interactions found for the provided genes.")
}

# --- 3. Process Graph Metadata ---
edges <- interactions %>%
  select(preferredName_A, preferredName_B, score) %>%
  rename(from = preferredName_A, to = preferredName_B)

g <- graph_from_data_frame(edges, directed = FALSE, vertices = all_genes_meta)

comp <- components(g)
g_sub <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])
cat("Filtered to giant component:", vcount(g_sub), "nodes,", ecount(g_sub), "edges.\n")

# --- 4. Visualize ---
cat("Generating visualization...\n")
set.seed(42)

plot_ppi <- ggraph(g_sub, layout = "stress") +
  geom_edge_link(aes(alpha = score), color = "grey70", show.legend = FALSE) +
  geom_node_point(aes(size = log10P, fill = type), shape = 21, color = "black", stroke = 0.6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold", color = "black", bg.color = "white", bg.r = 0.15) +
  scale_fill_manual(values = c("Lead SNP Gene" = "#ff7f00", "Pathway Driver" = "#1f78b4"), name = "Gene Type") +
  scale_size_continuous(range = c(4, 12), name = expression(-log[10](P[GWAS]))) +
  scale_alpha_continuous(range = c(0.1, 0.6)) +
  labs(
    title = "Systemic Lupus Erythematosus: Protein Interaction Network",
    subtitle = "57 Lead Loci Genes (STRING-DB Confidence > 0.400)",
    caption = "Node size reflects discovery significance. Edge transparency reflects interaction confidence."
  ) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right")

ggsave(OUTPUT_PNG, plot_ppi, width = 10, height = 10, dpi = 300, bg = "white")
cat("PPI Network Update Complete. Files saved to:", FIG_DIR, "\n")
