# AGENTS.md — Project Memory & Constraints

## 🛑 The 10 Critical Rules (DO NOT BREAK)
1.  **Lambda GC = 0.92**: The validated genomic inflation factor is 0.917. Never attempt to "correct" or "re-run" to force lambda = 1.0; this has been statistically justified for publication.
2.  **High-Confidence Criteria**: A locus is ONLY "High-Confidence" if it meets BOTH: **LAVA rho > 0.30** AND **COLOC PP4 > 0.80**.
3.  **Spanish Replication**: The final validation gate requires **P < 0.05** and **concordant direction** in the Spanish cohort.
4.  **Performance First**: Always use R `data.table` and `vroom` for summary statistics. Avoid `read.table` or base R data frames for large GWAS files.
5.  **No "Ghost" Software**: Do NOT reference METAL, SuSiE, or LDSC in the manuscript. All analysis is native R-based (IVW meta-analysis, coloc.abf, LAVA).
6.  **Reference Panel**: Always use **1000G Phase 3 European (EUR)** as the LD reference.
7.  **RSID Integrity**: `master_results_table.tsv` must have zero duplicate RSIDs and zero NA P-values.
8.  **Colocalization Tissues**: Primary tissues are **GTEx v10 Whole Blood and Spleen**.
9.  **Manuscript Sync**: Every figure in `SLE_Publication_Package/` must match the statistics in `results/master_results_table.tsv`.
10. **Aesthetic Standard**: All publication figures must use the established "Premium Minimalist" style (Black/White or subtle HSL palettes).

## 📂 File Map (Core Pipeline)
- `step1_meta_discovery.R`: Custom IVW fixed-effects meta-analysis (Bentham + FinnGen).
- `step2_spanish_replication.R`: Independent North-to-South European validation.
- `step5_figures.R`: Core GWAS visualizations (Manhattan, QQ).
- `step10_fgsea_ranked.R`: Simes-aggregated pathway enrichment (Reactome/Hallmark).
- `step20_precise_novelty.R`: Multi-database novelty audit (GWAS Catalog + PubMed).
- `step21_consolidate_final_results.R`: The "Source of Truth" table generator.
- `step22_lava_analysis.R`: Local genetic correlation analysis.
- `step24_coloc_analysis.R`: Bayesian colocalization (coloc.abf).
- `step26_phewas_lookup.R`: GWAS Catalog API-based pleiotropy audit.
- `step27_pleiotropy_plot.R`: High-aesthetic linear pleiotropy map.

## 🛠️ Common Task Recipes
- **Update Master Results**: `Rscript scripts/step21_consolidate_final_results.R`
- **Regenerate Publication Figures**: `Rscript scripts/step5_figures.R` followed by the `create_..._publication.R` specific scripts.
- **Run New Colocalization**: Update `step24_coloc_analysis.R` with the new target gene and run.
- **Check for Novelty**: `Rscript scripts/step20_precise_novelty.R`

## 📊 Environment Table
| Variable | Value | Description |
| :--- | :--- | :--- |
| `P_DISCOVERY` | 5e-8 | Genome-wide significance threshold |
| `P_REPLICATION`| 0.05 | Nominal significance for validation |
| `LAVA_RHO` | 0.30 | Threshold for high-confidence local rG |
| `COLOC_PP4` | 0.80 | Threshold for high-confidence colocalization |
| `LD_REF` | 1000G_P3_EUR | Reference panel for LD calculations |

## 📜 Change Log (Recent Highlights)
- **2026-04-25**: Finalized "Minimalist B&W" Pipeline Figure for supplementary.
- **2026-04-25**: Corrected Manuscript Methods (Removed METAL/LDSC/MAGMA references).
- **2026-04-24**: Resolved all RSID duplicates and NA P-values in results.
- **2026-04-23**: Finalized the 46-locus architecture.
