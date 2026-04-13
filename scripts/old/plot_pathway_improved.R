#!/usr/bin/env Rscript
# plot_pathway_improved.R
# Aesthetic improvement of fgsea pathway results for publication.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

# --- Configuration ---
DIR_IN  <- "results/pathway"
DIR_OUT <- "figures"
FILE_IN <- file.path(DIR_IN, "fgsea_significant_all.tsv")
OUTPUT_FIG <- file.path(DIR_OUT, "Fig_2_pathway_fgsea_v2.png")

# --- Load and Process Data ---
cat("Loading data from", FILE_IN, "...\n")
dt <- fread(FILE_IN)

# Function to clean pathway names
clean_pathway <- function(x) {
  x <- gsub("^(GOBP|REACTOME|HALLMARK|GSE\\d+)_", "", x)
  x <- gsub("_", " ", x)
  x <- str_to_sentence(x)
  # Capitalize specific acronyms
  x <- gsub("Mhc", "MHC", x)
  x <- gsub("Sle", "SLE", x)
  x <- gsub("Nfkb", "NF-kB", x)
  x <- gsub("Dc dn", "DC (DN)", x)
  return(x)
}

# Function to select top driver genes
get_drivers <- function(edge_str) {
  genes <- strsplit(edge_str, ";")[[1]]
  if (length(genes) > 3) {
    return(paste0(paste(genes[1:3], collapse=", "), "..."))
  }
  return(paste(genes, collapse=", "))
}

# Process for plotting
plot_dt <- dt %>%
  group_by(database) %>%
  # Select top 6 per database by significance
  arrange(padj) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  mutate(
    clean_name = sapply(pathway, clean_pathway),
    drivers = sapply(leadingEdge, get_drivers),
    # Combine name and drivers for a rich label
    rich_label = paste0(clean_name, "\n(", drivers, ")"),
    log10padj = -log10(padj),
    # Truncate label if too long
    rich_label = str_trunc(rich_label, 100)
  )

# Fix database names for faceting
plot_dt$database <- case_when(
  plot_dt$database == "GO_BP" ~ "GO Biological Process",
  plot_dt$database == "Hallmark" ~ "MSigDB Hallmark",
  TRUE ~ plot_dt$database
)

# Set factor order by NES within database
plot_dt <- plot_dt %>%
  arrange(database, NES) %>%
  mutate(rich_label = factor(rich_label, levels = rich_label))

# --- Generate Plot ---
cat("Generating improved bubble plot...\n")

p <- ggplot(plot_dt, aes(x = NES, y = rich_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # Add segments for visual link
  geom_segment(aes(x = 0, xend = NES, y = rich_label, yend = rich_label), 
               color = "grey90", linewidth = 0.5) +
  geom_point(aes(size = log10padj, fill = NES), shape = 21, color = "white", stroke = 0.5) +
  scale_size_continuous(range = c(5, 12), name = expression(-log[10](P[adj]))) +
  # Use a premium diverging palette
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", 
                       midpoint = 0, name = "NES") +
  facet_wrap(~database, scales = "free_y", ncol = 1) +
  labs(title = "Pathway Enrichment: SLE European Meta-Analysis",
       subtitle = "Pre-ranked fgsea on gene-level p-values (Simes method) | Top discoveries per database",
       x = "Normalised Enrichment Score (NES)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey95"),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey98", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9, lineheight = 0.8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey40", size = 10),
    panel.spacing = unit(1, "lines")
  )

# Save high-res
cat("Saving to", OUTPUT_FIG, "...\n")
# Higher height because of 4 facets
ggsave(OUTPUT_FIG, p, width = 10, height = 12, dpi = 300, bg = "white")

cat("Done.\n")
