#!/usr/bin/env Rscript
# Top loci P-value plot with gene and SNP labels

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

loci <- fread("results/top_loci_summary_table.tsv")

# Top 10 loci for plotting
top_loci <- loci %>%
  mutate(Gene_clean = sapply(strsplit(Gene, " / "), function(x) tail(x, 1)),
         label = paste0(Gene_clean, "\n(", RSID, ")")) %>%
  arrange(P_meta) %>%
  head(10)

# Plot with gene/SNP labels
p <- ggplot(top_loci, aes(x = reorder(label, -log10(P_meta)), y = -log10(P_meta))) +
  geom_point(aes(size = abs(BETA), color = Replicated), alpha = 0.8, shape = 16) +
  scale_color_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#9E9E9E"), name = "Replicated") +
  scale_size_continuous(name = "Effect Size |BETA|", range = c(4, 10)) +
  coord_flip() +
  labs(title = "Top 10 SLE-Associated Loci",
       subtitle = "Gene (rsID) with discovery meta -log10(P)",
       x = "", y = expression(-log[10](P[meta]))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10, face = "italic", color = "black"),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/top_loci_labeled.png", p, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("figures/top_loci_labeled.pdf", p, width = 10, height = 8)

cat("Top loci plot saved to figures/top_loci_labeled.png/pdf\n")