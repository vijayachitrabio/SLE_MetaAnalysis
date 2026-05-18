# Reproducibility steps

Run commands from the repository root unless otherwise noted.

## 1. Environment

Install the R and Python dependencies listed in `requirements.txt` and the package imports at the top of each R script. A Dockerfile is provided for a containerized environment.

## 2. Core analysis

The ordered analysis scripts are in `scripts/`:

```bash
Rscript scripts/step1_meta_discovery.R
Rscript scripts/step2_spanish_replication.R
Rscript scripts/step3_qc_summary_tables.R
Rscript scripts/step4_annotation.R
Rscript scripts/step5_figures.R
Rscript scripts/step6_pathways_eqtl.R
Rscript scripts/step8_sensitivity.R
Rscript scripts/step9_eqtls.R
Rscript scripts/step10_fgsea_ranked.R
Rscript scripts/step11_locus_plots.R
Rscript scripts/step12_novelty_check.R
Rscript scripts/step13_eqtl_plots.R
Rscript scripts/step16_refined_locus_analysis.R
Rscript scripts/step17_proper_ld_clumping.R
Rscript scripts/step18_ld_based_locus.R
Rscript scripts/step19_ld_confirmation.R
Rscript scripts/step20_precise_novelty.R
Rscript scripts/step21_consolidate_final_results.R
Rscript scripts/step22_lava_analysis.R
Rscript scripts/step23_consolidate_lava.R
Rscript scripts/step24_coloc_analysis.R
Rscript scripts/step25_final_synthesis.R
Rscript scripts/step26_phewas_lookup.R
Rscript scripts/step27_pleiotropy_plot.R
Rscript scripts/step28_susie_finemapping.R
Rscript scripts/step29_statistical_validation.R
```

## 3. Publication figure formatting

After the core result tables are present under `results/`, regenerate the final edited manuscript figures with:

```bash
Rscript scripts/publication_figures/make_figure_2_lava_prioritization.R
Rscript scripts/publication_figures/make_figure_4_replication_concordance.R
Rscript scripts/publication_figures/make_figure_5A_pathway_enrichment.R
```

These scripts write to `manuscript/figures/main/`.

## 4. Manuscript outputs

Curated manuscript outputs are organized as:

- `manuscript/figures/main/`: Figures 1-7, including Figure 5A and Figure 5B.
- `manuscript/figures/supplementary/`: Supplementary Figures 1-5.
- `manuscript/tables/main/`: Tables 1-4.
- `manuscript/tables/supplementary/`: Supplementary Tables 1-9.

The authoritative figure and table inventory is `list.txt`.
