# Deployment & Operations Guide

## 💻 Local Setup (macOS / Linux)
1. **R Environment**: Ensure R version >= 4.0.0 is installed.
2. **Library Installation**:
   ```R
   install.packages(c("data.table", "vroom", "dplyr", "ggplot2", "cowplot", "httr", "jsonlite"))
   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install(c("coloc", "fgsea"))
   ```
3. **Python Environment**:
   ```bash
   pip install pandas matplotlib numpy pathlib
   ```

## 📦 Data Prerequisites
The pipeline expects the following raw files (not included in the repo due to size):
- `/data/bentham_2015_gwas.tsv`
- `/data/finngen_r12_sle.tsv`
- `/data/spanish_replication_summary.tsv`
- `/reference/1000G_P3_EUR/` (Standard PLINK files)

## 🚀 Execution Flow
To run the full pipeline from scratch:
1.  **Run Discovery**: `Rscript scripts/step1_meta_discovery.R`
2.  **Annotate**: `Rscript scripts/step4_annotation.R`
3.  **Run LAVA**: `Rscript scripts/step22_lava_analysis.R`
4.  **Consolidate**: `Rscript scripts/step21_consolidate_final_results.R`

## 🛠️ Adding a New Page/Analysis
1. Create a new script in `scripts/` following the `stepXX_name.R` naming convention.
2. Update the File Map in `CLAUDE.md`.
3. If the script generates a result, ensure it appends or updates `results/master_results_table.tsv`.

## ⚠️ Troubleshooting
- **Memory Errors**: GWAS files are large. If R crashes, ensure you are using `fread` from `data.table` and consider running on a machine with >= 32GB RAM.
- **API Limits**: `step26_phewas_lookup.R` queries the GWAS Catalog. If you see 429 errors, increase the `Sys.sleep(0.5)` value.
- **Lambda Discrepancies**: If your lambda changes significantly, verify that you are calculating it on the full summary statistics *before* any p-value filtering.
