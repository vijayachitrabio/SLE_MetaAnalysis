# Manuscript Figures and Tables - SLE Meta-Analysis

## MAIN MANUSCRIPT FIGURES

| Figure | Description | Panel | Files to Combine |
|--------|-------------|-------|------------------|
| **Fig 1** | Meta-analysis Overview | A: Pipeline flowchart, B: Study characteristics | project_pipeline.dot → render |
| **Fig 2** | GWAS Results | A: Manhattan plot, B: QQ plot | manhattan_plot.png + qq_plot.png |
| **Fig 3** | Top Loci Forest Plots | Top 15 lead SNPs | forest_plot_top15.png |
| **Fig 4** | PPI Network | 23 nodes, 65 edges | Fig_3_PPI_network.png + .pdf |
| **Fig 5** | Pathway Enrichment | A: Dot plot, B: Bar plot | pathway_enrichment_combined.pdf |
| **Fig 6** | eQTL Multi-tissue | A: Dot heatmap, B: Clustered heatmap | Fig_4_eQTL_heatmap.png + Fig_4_eQTL_clustered_heatmap.png |
| **Fig 7** | Therapeutic Mapping | Drug repurposing candidates | therapeutic_mapping.pdf |

## MAIN MANUSCRIPT TABLES

| Table | Description | Source |
|-------|-------------|--------|
| **Table 1** | Study cohorts and sample sizes | sample_size_summary.tsv |
| **Table 2** | QC metrics (λ, inflation) | lambda_qc_report.tsv + inflation_summary.tsv |
| **Table 3** | Top loci with annotation | top_loci_summary_table.tsv |
| **Table 4** | Replication results | sensitivity/SA3_replication_power.tsv |

---

## SUPPLEMENTARY FIGURES (S-Figures)

| Figure | Description | Files to Combine |
|--------|-------------|------------------|
| **S-Fig 1** | Locus zoom plots (5 loci) | locus_plots/*.png (IRF5, STAT4, ITGAM, TNIP1, HLA-MHC) |
| **S-Fig 2** | Heterogeneity by locus | SA1_heterogeneity.png + .pdf |
| **S-Fig 3** | HLA pairwise LD | SA2_hla_distance.png |
| **S-Fig 4** | Replication power analysis | SA3_replication.png + .pdf |
| **S-Fig 5** | Novelty classification | novelty_classification.png + .pdf |
| **S-Fig 6** | Top 30 loci bar plot | top_loci_summary.png + .pdf |
| **S-Fig 7** | Pathway NES bubbles | pathway_nes_bubbles.png |
| **S-Fig 8** | Pathway bar plot | pathway_enrichment_barplot.png |

---

## SUPPLEMENTARY TABLES (S-Tables)

| Table | Description | Source |
|-------|-------------|--------|
| **S-Table 1** | Full discovery meta-analysis results | discovery_meta_results.tsv (full) |
| **S-Table 2** | Gene-level P-values | gene_level_pvalues.tsv |
| **S-Table 3** | Full pathway enrichment (all terms) | pathway_enrichment_results.tsv |
| **S-Table 4** | Filtered significant pathways (top 25) | pathway_enrichment_filtered.tsv |
| **S-Table 5** | eQTL summary by tissue | eqtl_summary.tsv |
| **S-Table 6** | Heterogeneity (I², τ²) for lead SNPs | heterogeneity_lead_snps.tsv |
| **S-Table 7** | SA1: Random vs Fixed effects | sensitivity/SA1_random_vs_fixed.tsv |
| **S-Table 8** | SA2: HLA pairwise distances | sensitivity/SA2_HLA_pairwise.tsv |
| **S-Table 9** | SA3: Replication power | sensitivity/SA3_replication_power.tsv |
| **S-Table 10** | Therapeutic drug status summary | drug_status_summary.tsv |
| **S-Table 11** | Novel loci classification | Based on Known_SLE in top_loci_summary_table.tsv |

---

## PANEL COMBINATIONS (for single figures)

### Fig 1: Overview
- A: Pipeline (project_pipeline.dot rendered)
- B: Sample size summary (from sample_size_summary.tsv)

### Fig 2: GWAS Results  
- A: manhattan_plot.png
- B: qq_plot.png

### Fig 3: Top Loci
- Combine: forest_plot_top15.png

### Fig 4: PPI Network
- Single panel: Fig_3_PPI_network.png

### Fig 5: Pathway Enrichment
- A: pathway_enrichment.png (dot plot)
- B: pathway_enrichment_barplot.png (bar plot)
- Or: pathway_enrichment_combined.pdf

### Fig 6: eQTL
- A: Fig_4_eQTL_heatmap.png (dot plot)
- B: Fig_4_eQTL_clustered_heatmap.png (clustered heatmap)

### Fig 7: Therapeutic
- Single panel: therapeutic_mapping.png

---

## Notes
- All PDFs are publication-ready vector format
- PNGs at 300 DPI suitable for print
- Locus plots available for: IRF5, STAT4, ITGAM, TNIP1, HLA-MHC region