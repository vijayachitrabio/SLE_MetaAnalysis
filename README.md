# SLE Meta-Analysis: European Ancestry and Functional Integration

This repository contains the analytical pipeline and scripts for a fixed-effects meta-analysis of Systemic Lupus Erythematosus (SLE) GWAS in European ancestry populations, with downstream functional genomics integration.

## Project Overview

We performed an inverse-variance weighted (IVW) meta-analysis of two European GWAS cohorts (**Bentham 2015** and **Julià 2018**), totaling 13,352 individuals (4,835 cases, 8,517 controls). We identified **13 independent genome-wide significant loci** ($P < 5 \times 10^{-8}$).

### Key Features
- **Meta-Analysis**: Robust IVW fixed-effects model with per-cohort genomic control.
- **Independent Pruning**: Conservative physical distance pruning (1MB for MHC, 500kb elsewhere).
- **Replication**: Validation in **FinnGen Release 12** (L12_LUPUS).
- **Functional Integration**:
    - Multi-tissue cis-eQTL mapping using **GTEx v10**.
    - Immune cell eQTL validation using the **DICE** database (B cells, T cells, Monocytes).


## Repository Structure

- `scripts/`: R scripts for meta-analysis, QC, pathway enrichment, and functional integration.
- `results/`: Key summary tables (lead SNPs, functional evidence, annotation).

## Usage

Scripts are intended to be run sequentially:
1. `scripts/sle_meta_analysis.R`: Core meta-analysis and global QC.
2. `scripts/step2_annotation_regional.R`: SNP-to-gene mapping (proximity and eQTL).
3. `scripts/sle_ld_clumping_v2.R`: Identification of independent lead SNPs.
4. `scripts/step4_pathway_publication.R`: Gene-level aggregation and pathway enrichment (GSEA/gProfiler).
5. `scripts/get_gtex_eqtls.R`: Integration with GTEx API for top loci.
6. `scripts/step6_coloc_summary.R`: Functional evidence consolidation.

## Functional Evidence Highlights

Our analysis identified five loci with **Strong** multi-source functional evidence (GTEx + DICE), directly implicating the following genes in SLE pathogenesis:
- **BTN3A2**
- **HLA-DQA2 / DRB6**
- **IRF5**
- **FAM167A / BLK**
- **C4A**

## Data Sources & Availability

**Discovery Meta-Analysis (European Ancestry):**
- **Bentham et al. (2015)**: *Genetic association study of systemic lupus erythematosus.* [GCST003156](https://www.ebi.ac.uk/gwas/studies/GCST003156)
- **Julià et al. (2018)**: *Genome-wide association study identifies novel SLE risk loci.*
**Independent Replication (European Ancestry):**
- **FinnGen Release 12**: Systemic Lupus Erythematosus (Phenocode: `L12_LUPUS`).
    - *Cases: ~850 | Controls: ~418,000*
    - [FinnGen R12 Phenotype Page](https://r12.finngen.fi/pheno/L12_LUPUS)


---
*Maintained by [vijayachitrabio](https://github.com/vijayachitrabio)*
