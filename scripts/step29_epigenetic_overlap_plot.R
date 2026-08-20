#!/usr/bin/env Rscript
# scripts/step29_epigenetic_overlap_plot.R
# Generate Premium Minimalist plot for Epigenetic Overlap Results

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

cat("Loading overlap results...\n")
res <- fread("results_extracted/exploratory_epigenetic_overlap.tsv")

# Filter out MHC
res <- res[Is_MHC == FALSE]
total_loci <- nrow(res)

# Calculate percentages for each cell type + aggregate
calc_stats <- function(df, col_exact, col_window) {
    pct_exact <- (sum(df[[col_exact]]) / total_loci) * 100
    pct_window <- (sum(df[[col_window]]) / total_loci) * 100
    return(c(pct_exact, pct_window))
}

stats_pbmc <- calc_stats(res, "E062_PBMC_Exact", "E062_PBMC_Window10kb")
stats_bcell <- calc_stats(res, "E032_B_cell_Exact", "E032_B_cell_Window10kb")
stats_tcell <- calc_stats(res, "E034_T_cell_Exact", "E034_T_cell_Window10kb")
stats_mono <- calc_stats(res, "E029_Monocyte_Exact", "E029_Monocyte_Window10kb")
stats_any <- calc_stats(res, "Any_Immune_Exact", "Any_Immune_Window10kb")

plot_data <- data.frame(
    Cell_Type = factor(c("PBMC (E062)", "B Cell (E032)", "T Cell (E034)", "Monocyte (E029)", "Any Immune Cell"),
                       levels = c("PBMC (E062)", "B Cell (E032)", "T Cell (E034)", "Monocyte (E029)", "Any Immune Cell")),
    Exact = c(stats_pbmc[1], stats_bcell[1], stats_tcell[1], stats_mono[1], stats_any[1]),
    Window = c(stats_pbmc[2], stats_bcell[2], stats_tcell[2], stats_mono[2], stats_any[2])
)

plot_data_long <- pivot_longer(plot_data, cols = c("Exact", "Window"), names_to = "Overlap_Type", values_to = "Percentage")
plot_data_long$Overlap_Type <- factor(plot_data_long$Overlap_Type, levels = c("Exact", "Window"), 
                                      labels = c("Exact Lead-SNP", "Proximal Window (±10kb)"))

# Premium Minimalist Aesthetics (per AGENTS.md rules)
cat("Generating plot...\n")
p <- ggplot(plot_data_long, aes(x = Cell_Type, y = Percentage, fill = Overlap_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Exact Lead-SNP" = "#404040", "Proximal Window (±10kb)" = "#A0A0A0")) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(
        x = "Immune Cell Annotation (Roadmap Epigenomics)",
        y = "Percentage of Non-MHC Loci (%)",
        fill = "Overlap Criterion",
        title = "SLE Loci Overlap with Immune Enhancer/Promoter States"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 15)),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black", linewidth = 0.5),
        plot.margin = margin(20, 20, 20, 20)
    )

# Save to SLE_Publication_Package
dir.create("SLE_Publication_Package", showWarnings = FALSE)
output_pdf <- "SLE_Publication_Package/Supplementary_Figure_Epigenetic.pdf"
output_png <- "SLE_Publication_Package/Supplementary_Figure_Epigenetic.png"

ggsave(output_pdf, plot = p, width = 8, height = 6, dpi = 300)
ggsave(output_png, plot = p, width = 8, height = 6, dpi = 300)

cat(sprintf("Success! Saved figure to %s and %s\n", output_pdf, output_png))
