# Technical Architecture: The SLE Meta-Analysis Pipeline

## 🏗️ System Overview
The pipeline is a modular, R-based framework designed for high-throughput GWAS meta-analysis and causal prioritization. It consists of 29 logic steps, divided into four primary phases.

## 🧵 Phase 1: Discovery & Quality Control
1.  **Harmonization**: Alignment of effect alleles across Bentham and FinnGen.
2.  **Meta-Analysis**: Native R implementation of Inverse Variance Weighted (IVW) fixed-effects.
3.  **GC Correction**: Application of lambda-based correction ($\lambda = 0.917$).
4.  **Clumping**: PLINK-style LD clumping using 1000G EUR reference ($r^2 < 0.1$, distance = 500kb).

## 🧵 Phase 2: Independent Validation
1.  **Spanish Replication**: Look-up of 47 lead SNPs in the Southern European cohort.
2.  **Sign Concordance**: Filtering for SNPs with matching beta signs in discovery and replication.
3.  **Heterogeneity**: Calculation of Cochran’s Q and $I^2$ to assess cohort consistency.

## 🧵 Phase 3: Causal Prioritization (LAVA & COLOC)
1.  **LAVA (Local Genetic Correlation)**: Assessing whether genetic signals in Whole Blood/Spleen align with the SLE signal.
2.  **Bayesian Colocalization (coloc.abf)**: Quantifying the probability of a shared causal variant (PP4).
3.  **Tiering**: Categorization into High-Confidence (LAVA rho > 0.3, PP4 > 0.8) vs. Evidence-Supported.

## 🧵 Phase 4: Interpretation & Translation
1.  **FGSEA**: Gene-set enrichment against MSigDB collections.
2.  **PheWAS**: API-based pleiotropy mapping against the GWAS Catalog.
3.  **Therapeutic Mapping**: Prioritization of drugs based on FDA approval and clinical trial status.

## 🌳 Annotated File Tree
```text
.
├── scripts/
│   ├── step1_meta_discovery.R      # Primary meta-analysis logic
│   ├── step22_lava_analysis.R     # Local genetic correlation (LAVA)
│   ├── step24_coloc_analysis.R    # Bayesian colocalization (COLOC)
│   └── create_..._publication.R   # High-aesthetic figure generation
├── results/
│   ├── master_results_table.tsv   # THE SOURCE OF TRUTH (All loci)
│   └── lava_results_final.tsv     # Processed LAVA correlations
└── SLE_Publication_Package/       # FINAL MANUSCRIPT ASSETS
    ├── Tables/                    # Formal TSV/CSV tables
    └── Figures/                   # High-res PNG/PDF assets
```

## 💻 Tech Stack
- **Languages**: R (90%), Python (10% for complex viz).
- **Core Libraries (R)**: `data.table`, `vroom`, `ggplot2`, `coloc`, `LAVA`, `fgsea`.
- **Core Libraries (Python)**: `pandas`, `matplotlib`, `pathlib`.
- **Reference Panel**: 1000 Genomes Phase 3 (EUR subset).
