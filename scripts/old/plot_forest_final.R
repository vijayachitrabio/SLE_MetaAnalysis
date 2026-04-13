#!/usr/bin/env Rscript
# plot_forest_final.R
# Standalone script to generate a publication-quality Forest plot
# with replication symbols, diamond sizing by P-value, and status coloring.

library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

# --- Configuration ---
DIR_RES <- "results"
INPUT_FILE <- file.path(DIR_RES, "european_replication_concordance_table.tsv")
OUTPUT_PNG <- file.path(DIR_RES, "final_forest_plot_publication.png")
OUTPUT_PDF <- file.path(DIR_RES, "final_forest_plot_publication.pdf")

# --- Load and Prepare Data ---
cat("Loading data from", INPUT_FILE, "...\n")
dt <- fread(INPUT_FILE)

# Prepare columns for plotting
plot_df <- dt %>%
  mutate(
    # Odds Ratio and CI
    OR = exp(BETA_meta),
    lower_CI = exp(BETA_meta - 1.96 * SE_meta),
    upper_CI = exp(BETA_meta + 1.96 * SE_meta),
    
    # Gene and Variant labels
    clean_gene = ifelse(Gene == "" | is.na(Gene), SNP, Gene),
    variant_label = paste0(clean_gene, " (", SNP, ")"),
    effect_allele = paste0(A1, ">", A2),
    label_full = paste0(variant_label, "\n", effect_allele),
    chr_pos = paste0("Chr", CHR, ":", format(BP, big.mark=",")),
    
    # Sizing
    log10p = -log10(P_meta),
    
    # Replication Status and Symbols
    status = case_when(
      replication_status == "palindromic_excluded" ~ "Palindromic excluded",
      replicated == TRUE ~ "Replicated",
      replicated == FALSE ~ "Not replicated",
      TRUE ~ "Not tested"
    ),
    status_symbol = case_when(
      replication_status == "palindromic_excluded" ~ "P",
      replicated == TRUE ~ "\u2713", # Check mark
      replicated == FALSE ~ "\u2717", # X mark
      TRUE ~ ""
    )
  )

# Order by genomic position (Chromosome and BP)
plot_df <- plot_df %>%
  arrange(CHR, BP) %>%
  mutate(label_full = factor(label_full, levels = rev(label_full)))

# Factor for colors
plot_df$status <- factor(plot_df$status, levels = c("Replicated", "Not replicated", "Palindromic excluded"))

# --- Generate Plot ---
cat("Generating publication forest plot...\n")

# Main plot
p_main <- ggplot(plot_df, aes(x = OR, y = label_full)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2, color = "grey40", linewidth = 0.5) +
  geom_point(aes(size = log10p, color = status), shape = 18) +
  # Using trans = "log" but NO limits here to avoid removing points at the edge
  scale_x_continuous(trans = "log", breaks = c(0.5, 0.7, 1, 1.3, 1.5, 1.8, 2.0),
                     name = "Odds Ratio (95% CI)") +
  scale_color_manual(values = c("Replicated" = "#2166ac", 
                                "Not replicated" = "#d73027", 
                                "Palindromic excluded" = "grey60")) +
  scale_size_continuous(range = c(3, 11), name = expression(-log[10](P[meta]))) +
  labs(title = "SLE European Meta-Analysis: 13 Lead Variants",
       subtitle = "Discovery: Bentham + Julià (4,835 cases)  |  Replication: FinnGen R12 (850 cases)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey95"),
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(10, 80, 10, 10) # Add large right margin for symbols
  )

# Add the symbol column on the right
# We use x = 2.4 to place it outside the main plot area
p_final <- p_main +
  geom_text(aes(x = 2.4, label = status_symbol, color = status), 
            hjust = 0.5, size = 6, fontface = "bold", show.legend = FALSE) +
  # Use coord_cartesian to set display limits WITHOUT removing data
  coord_cartesian(xlim = c(0.55, 2.1), clip = "off") + 
  annotate("text", x = 2.4, y = length(unique(plot_df$label_full)) + 0.8, 
           label = "FinnGen\nStatus", fontface = "bold", size = 3.5, vjust = 0)

# Save high-res
cat("Saving high-resolution plot to PNG and PDF formats...\n")
ggsave(OUTPUT_PNG, p_final, width = 11, height = 9, dpi = 300, bg = "white")
ggsave(OUTPUT_PDF, p_final, width = 11, height = 9, device = "pdf", bg = "white")

cat("Done.\n")
