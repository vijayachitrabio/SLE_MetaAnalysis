#!/usr/bin/env Rscript
# plot_gprofiler_improved.R
# aesthetic improvement of g:Profiler results for publication supplement.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

# --- Configuration ---
DIR_IN  <- "results/pathway"
DIR_OUT <- "figures"
FILE_IN <- file.path(DIR_IN, "gprofiler_results.tsv")
OUTPUT_FIG <- file.path(DIR_OUT, "Fig_S2_pathway_gprofiler_v2.png")

# --- Load and Process Data ---
cat("Loading data from", FILE_IN, "...\n")
dt <- fread(FILE_IN)

# Function to select top driver genes (intersections)
get_drivers <- function(intersection_str) {
  genes <- strsplit(intersection_str, ",")[[1]]
  if (length(genes) > 5) {
    return(paste0(paste(genes[1:5], collapse=", "), "..."))
  }
  return(paste(genes, collapse=", "))
}

# Process for plotting
plot_dt <- dt %>%
  filter(significant == TRUE) %>%
  group_by(source) %>%
  # Select top 10 per source by p-value
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    clean_name = str_to_sentence(term_name),
    drivers = sapply(intersection, get_drivers),
    # Rich label
    rich_label = paste0(clean_name, " (", term_id, ")\n", drivers),
    log10p = -log10(p_value)
  )

# Set factor order
plot_dt <- plot_dt %>%
  arrange(source, log10p) %>%
  mutate(rich_label = factor(rich_label, levels = unique(rich_label)))

# --- Generate Plot ---
cat("Generating improved dotplot...\n")

p <- ggplot(plot_dt, aes(x = log10p, y = rich_label)) +
  geom_segment(aes(x = 0, xend = log10p, y = rich_label, yend = rich_label), 
               color = "grey90") +
  geom_point(aes(size = intersection_size, color = precision)) +
  scale_size_continuous(range = c(3, 8), name = "Intersection Size") +
  scale_color_viridis_c(option = "magma", end = 0.9, name = "Precision") +
  facet_wrap(~source, scales = "free_y", ncol = 1) +
  labs(title = "Overrepresentation Analysis: g:Profiler",
       subtitle = "Significant GO and Reactome terms from lead SNP regions",
       x = expression(-log[10](p-value)),
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey98", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8, lineheight = 0.8),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# Save
cat("Saving to", OUTPUT_FIG, "...\n")
ggsave(OUTPUT_FIG, p, width = 11, height = 13, dpi = 300, bg = "white")

cat("Done.\n")
