#!/usr/bin/env Rscript
# scripts/plot_ppi_network.R
# Generate a high-quality Protein-Protein Interaction (PPI) network using lead SNP genes 
# and pathway driver genes, with dynamic data loading and premium aesthetics.

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

# --- Configuration ---
RESULTS_DIR <- "results"
FIG_DIR     <- "figures"
dir.create(FIG_DIR, showWarnings = FALSE)

FILE_ANNOT    <- file.path(RESULTS_DIR, "gene_annotation.tsv")
FILE_PATHWAY  <- file.path(RESULTS_DIR, "pathway/gprofiler_results.tsv")
OUTPUT_PNG    <- file.path(FIG_DIR, "Fig_3_PPI_network_improved.png")
OUTPUT_PDF    <- file.path(FIG_DIR, "Fig_3_PPI_network_improved.pdf")

cat("=== Building Enhanced PPI Network ===\n")

# --- 1. Load Data Dynamically ---
if (!file.exists(FILE_ANNOT)) stop("Gene annotation file not found.")
if (!file.exists(FILE_PATHWAY)) stop("Pathway results file not found.")

# A. Lead SNP Genes
annot_dt <- fread(FILE_ANNOT)
# Filter for actual gene names (remove "HLA-upstream" or empty/NA)
seed_genes_df <- annot_dt[Gene != "" & !is.na(Gene) & !grepl("upstream", Gene)]

# Resolve "TNF/HLA" style names and handle duplicates
seed_genes_df <- seed_genes_df %>%
  mutate(Gene = str_split(Gene, "/") %>% sapply(`[`, 1)) %>% # Take first gene if split by /
  group_by(Gene) %>%
  summarize(log10P = max(-log10(P), na.rm = TRUE), .groups = "drop")

# B. Pathway Driver Genes (Extract from top 10 significant pathways)
path_dt <- fread(FILE_PATHWAY)
driver_list <- path_dt %>%
  filter(significant == TRUE) %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  pull(intersection) %>%
  paste(collapse = ",") %>%
  str_split(",") %>%
  unlist() %>%
  unique()

# C. Merge and Label
all_genes_meta <- data.frame(Gene = unique(c(seed_genes_df$Gene, driver_list))) %>%
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
  required_score = 700, # High confidence
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
V(g)$degree <- degree(g)

# Only keep the main connected components if fragmented
comp <- components(g)
g_sub <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])
cat("Filtered to giant component:", vcount(g_sub), "nodes,", ecount(g_sub), "edges.\n")

# --- 4. Visualize with Premium Aesthetics ---
cat("Generating visualization...\n")

# Set seed for reproducible layout
set.seed(123)

plot_ppi <- ggraph(g_sub, layout = "stress") +
  # Edges with score-based alpha and color
  geom_edge_link(aes(alpha = score), color = "grey70", show.legend = FALSE) +
  geom_edge_arc(aes(alpha = score), color = "grey80", strength = 0.1, show.legend = FALSE) +
  
  # Nodes with GWAS P-value sizing and Type-based coloring
  geom_node_point(aes(size = log10P, fill = type), shape = 21, color = "black", stroke = 0.6) +
  
  # Reputable labels using geom_node_text
  geom_node_text(aes(label = name), 
                 repel = TRUE, size = 4, fontface = "bold",
                 color = "black", bg.color = "white", bg.r = 0.15) +
  
  # Premium color palette and sizing
  scale_fill_manual(values = c("Lead SNP Gene" = "#ff7f00", "Pathway Driver" = "#1f78b4"),
                    name = "Gene Type") +
  scale_size_continuous(range = c(5, 15), name = expression(-log[10](P[GWAS]))) +
  scale_alpha_continuous(range = c(0.1, 0.6)) +
  
  # Theme and Labels
  labs(
    title = "Systemic Lupus Erythematosus: Protein Interaction Network",
    subtitle = "Lead GWAS Loci and Enriched Pathway Drivers (STRING-DB High Confidence > 0.700)",
    caption = "Node size reflects GWAS discovery significance. Edge transparency reflects STRING interaction confidence."
  ) +
  theme_graph(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 20, margin = margin(b=5)),
    plot.subtitle = element_text(size = 12, color = "grey30", margin = margin(b=15)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save
cat("Saving improved figures...\n")
ggsave(OUTPUT_PNG, plot_ppi, width = 12, height = 11, dpi = 300, bg = "white")
ggsave(OUTPUT_PDF, plot_ppi, width = 12, height = 11, device = "pdf", bg = "white")

cat("PPI Network Update Complete. Files saved to:", FIG_DIR, "\n")
