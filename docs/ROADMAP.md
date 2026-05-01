# Roadmap & Manuscript Backlog

## 🟢 P0: Final Submission Readiness (This Week)
- [x] Finalize 47-locus master table (No duplicates/NAs).
- [x] Correct Manuscript Methods (Remove SuSiE/METAL/LDSC).
- [x] Generate B&W Supplementary Pipeline Figure.
- [ ] Run final sanity check on LAVA results against 1000G LD blocks.
- [ ] Finalize Table 1 (Cohort Demographics) in the manuscript.

## 🟡 P1: Supplemental Enhancements
- [ ] Generate individual LocusZoom-style plots for the 22 Novel loci.
- [ ] Create a "Spanish vs Discovery" forest plot for all 47 loci.
- [ ] Consolidate all script logs into a single `results/QC_audit_log.txt`.

## 🔵 P2: Future Extension Ideas
- **Polygenic Risk Score (PRS)**: Evaluate the performance of our 47-locus PRS in non-European populations.
- **Single-Cell eQTL**: Integrate the 15 high-confidence targets with immune single-cell data (e.g., OneK1K).
- **Mendelian Randomization (MR)**: Perform proteome-wide MR to identify causal circulating proteins.

## 🤖 AI Prompt Library
Use these prompts with Cursor/Claude to maintain the project:

### Prompt: Add a new locus
> "I have a new locus at CHR X, POS Y. Can you run the `step18_ld_based_locus.R` logic to define the 500kb window, pull the lead SNP, and check against `results/master_results_table.tsv` for existing overlaps?"

### Prompt: Style a publication figure
> "I need to visualize the results of [analysis] for publication. Please follow the `docs/DESIGN_SYSTEM.md` rules: use the `DejaVu Sans` font scale, avoid generic colors, and use the `cowplot` theme. Output should be 300 DPI PNG and PDF."

### Prompt: Sanity check results
> "Please scan `results/master_results_table.tsv`. Ensure there are exactly 47 unique RSIDs, zero NA p-values, and check if any P-values are exactly 0 (which should be converted to the smallest representable float)."
