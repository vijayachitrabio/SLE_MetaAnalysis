#!/usr/bin/env Rscript
# scripts/step7_ppi_network.R
# Generate a clean Protein-Protein Interaction (PPI) network

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
OUTPUT_PNG    <- file.path(FIG_DIR, "Fig_3_PPI_network.png")
OUTPUT_PDF    <- file.path(FIG_DIR, "Fig_3_PPI_network.pdf")

cat("=== Building PPI Network ===\n")

if (!file.exists(FILE_ANNOT)) stop("Gene annotation file not found.")

# Lead SNP Genes
annot_dt <- fread(FILE_ANNOT)
annot_dt <- annot_dt[Gene != "TBD" & Gene != "Intergenic" & !is.na(Gene)]

genes_df <- annot_dt %>%
  mutate(Gene = str_split(Gene, " / ") %>% sapply(function(x) x[length(x)])) %>% 
  group_by(Gene) %>%
  summarize(log10P = max(-log10(P_disco), na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "GWAS Gene")

cat("Genes for PPI:", nrow(genes_df), "\n")

# Fetch Interactions from STRING API
cat("Fetching PPI data from STRING API...\n")
string_api_url <- "https://string-db.org/api/json/network"
res <- GET(string_api_url, query = list(
  identifiers = paste(genes_df$Gene, collapse = "\r"),
  species = 9606,
  required_score = 400,
  caller_identity = "sle_meta_analysis"
))

if (status_code(res) != 200) {
  stop("Failed to fetch from STRING API. Status: ", status_code(res))
}

interactions <- fromJSON(content(res, as = "text"))
if (length(interactions) == 0) {
  stop("No interactions found.")
}

edges <- interactions %>%
  select(preferredName_A, preferredName_B, score) %>%
  rename(from = preferredName_A, to = preferredName_B)

g <- graph_from_data_frame(edges, directed = FALSE, vertices = genes_df)

# Keep largest connected component
comp <- components(g)
g_sub <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])
cat("Network:", vcount(g_sub), "nodes,", ecount(g_sub), "edges\n")

# Clean visualization
set.seed(42)
ppi <- ggraph(g_sub, layout = "fr") +
  geom_edge_link(aes(alpha = score), color = "grey60", show.legend = FALSE) +
  geom_node_point(aes(size = log10P, fill = type), shape = 21, color = "black", stroke = 0.5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "italic", 
                 color = "black", bg.color = "white", bg.r = 0.1) +
  scale_fill_brewer(palette = "Set1", name = "") +
  scale_size_continuous(range = c(3, 9), name = "-log10(P)") +
  theme_graph(base_family = "sans", background = "white") +
  theme(legend.position = "bottom", legend.text = element_text(size = 10))

ggsave(OUTPUT_PNG, ppi, width = 10, height = 10, dpi = 300, bg = "white")
ggsave(OUTPUT_PDF, ppi, width = 10, height = 10)
cat("PPI Network saved to figures/\n")
