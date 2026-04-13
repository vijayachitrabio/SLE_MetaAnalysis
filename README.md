# SLE Meta-Analysis: North-to-South European Geo-Validation

This repository contains the complete analytical pipeline and manuscript-quality results for a genome-wide association (GWAS) meta-analysis of Systemic Lupus Erythematosus (SLE) in European-ancestry populations.

## 1. Project Overview

### High-Impact Study Design: "North-to-South" Validation
The core novelty of this work lies in its continental genetic validation strategy. We leverage high-powered datasets from Northern/North-Eastern Europe for discovery, and an independent Southern European cohort for validation to characterize the stability of SLE genetic architecture across the European continent.

### 2. Analytical Statistics

- **Stage 1: Discovery Phase (Northern Europe)**
    - **Bentham 2015 (UK/Northern)**: 4,036 cases | 6,959 controls
    - **FinnGen R12 (M13_SLE)**: 1,306 cases | 418,172 controls
    - **Total Discovery**: **5,342 cases** identified **57 Independent Loci**.

- **Stage 2: Validation Phase (Southern Europe)**
    - **Julià 2018 (Spanish Only)**: 907 cases | 1,558 controls
    - **Validation Rate**: **26 loci (58%)** successfully replicated ($P < 0.05$).

- **Novelty Verification**
    - **31 Novel Loci (54%)** identified when cross-referencing with the official EBI GWAS Catalog for SLE.

---

## 3. Pipeline Structure (13 Modern Steps)

The repository is organized as a modular 13-step R pipeline:

| Step | Script | Key Output |
| :--- | :--- | :--- |
| **1-3** | [Scripts 1-3](scripts/) | Discovery Meta-analysis, Spanish Validation, and Results Synthesis. |
| **4** | `step4_annotation.R` | Mapping 57 loci to target protein-coding genes (GRCh38). |
| **5** | `step5_figures.R` | Manhattan, Q-Q, and Forest plots. |
| **6-7** | `step6_pathways_eqtl.R`, `step7_ppi_network.R` | gProfiler2 enrichment and STRING-DB interaction networks. |
| **8** | `step8_sensitivity.R` | Random vs Fixed effects, Replication Power, and HLA spacing. |
| **9** | `step9_eqtls.R` | Immune-tissue eQTL mapping via GTEx API. |
| **10** | `step10_fgsea_ranked.R` | Advanced ranked-list pathway analysis (NES). |
| **11** | `step11_locus_plots.R` | Regional Manhattan plots for the Top 5 loci (*STAT4*, *IRF5*, *HLA*, *ITGAM*, *TNIP1*). |
| **12** | `step12_novelty_check.R` | Formal "Novel vs Known" classification using `gwasrapidd`. |
| **13** | `step13_eqtl_plots.R` | Multi-tissue eQTL dot plot heatmap. |

---

## 4. Key Visualizations

The pipeline generates publication-ready figures in the `figures/` directory:
- **Main GWAS**: `figures/Manhattan_SLE_Meta.png` and `figures/QQ_Plot.png`.
- **Function**: `figures/Fig_3_PPI_network_new.png` and `figures/Fig_4_eQTL_heatmap.png`.
- **Regional**: `figures/locus_plots/` containing zoom-in association views for lead genes.

## 5. Getting Started

### Dependencies
Ensure the following R packages are installed:
`vroom`, `dplyr`, `data.table`, `biomaRt`, `fgsea`, `ggplot2`, `msigdbr`, `gwasrapidd`, `igraph`, `ggraph`.

### Execution
Scripts are intended to be run sequentially from the project root:
```bash
/usr/local/bin/Rscript scripts/step1_meta_discovery.R
...
/usr/local/bin/Rscript scripts/step13_eqtl_plots.R
```

---
*Maintained by [vijayachitramodhukur](https://github.com/vijayachitramodhukur)*
