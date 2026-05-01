library(dplyr)
library(readr)
library(stringr)
library(tidyr)

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

fmt_or_ci <- function(beta, se) {
  ifelse(
    is.na(beta) | is.na(se),
    "NA",
    sprintf(
      "%.2f (%.2f-%.2f)",
      exp(beta),
      exp(beta - 1.96 * se),
      exp(beta + 1.96 * se)
    )
  )
}

pick_gene <- function(effector_gene, gene) {
  effector_gene <- str_trim(replace_na(effector_gene, ""))
  gene <- str_trim(replace_na(gene, ""))
  ifelse(effector_gene != "", effector_gene, ifelse(gene != "", gene, "NA"))
}

master <- read_tsv("results/master_results_table.tsv", show_col_types = FALSE) %>%
  mutate(
    Top_Prioritized_Gene = pick_gene(Effector_Gene, Gene)
  ) %>%
  filter(Final_Assessment == "HIGH CONFIDENCE") %>%
  arrange(P_meta) %>%
  distinct(RSID, .keep_all = TRUE) %>%
  select(
    RSID, CHR, BP, P_meta, Novelty, Region, Top_Prioritized_Gene,
    LAVA_Status, Functional_Evidence, Best_PP4, Coloc_Tissue
  )

disco <- read_tsv("results/discovery_meta_results.tsv", show_col_types = FALSE) %>%
  select(RSID, OA, EA, BETA_meta, SE_meta, HetP, I2)

qc <- read_tsv("results/supplementary_qc_table.tsv", show_col_types = FALSE) %>%
  select(RSID, HetP_qc = HetP, I2_qc = I2)

therapy <- read_tsv("results/therapeutic_mapping_summary.tsv", show_col_types = FALSE) %>%
  mutate(
    therapy_text = case_when(
      !is.na(Therapeutic_Pathway) & Therapeutic_Pathway != "" &
        !is.na(Specific_Drug) & Specific_Drug != "" &
        !is.na(Drug_Status) & Drug_Status != "" ~
          paste0(
            "Pathway: ", Therapeutic_Pathway,
            "; candidate agents: ", Specific_Drug,
            " (", Drug_Status, ")"
          ),
      !is.na(Therapeutic_Pathway) & Therapeutic_Pathway != "" &
        !is.na(Drug_Class) & Drug_Class != "" ~
          paste0(
            "Pathway: ", Therapeutic_Pathway,
            "; drug class: ", Drug_Class
          ),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(therapy_text) & therapy_text != "") %>%
  group_by(RSID) %>%
  summarise(
    Therapy_Evidence = paste(unique(therapy_text), collapse = " | "),
    .groups = "drop"
  )

out <- master %>%
  left_join(disco, by = "RSID") %>%
  left_join(qc, by = "RSID") %>%
  left_join(therapy, by = "RSID") %>%
  mutate(
    HetP_final = coalesce(HetP_qc, HetP),
    I2_final = coalesce(I2_qc, I2),
    `chr:pos:A1:A2` = paste(CHR, BP, EA, OA, sep = ":"),
    `Lead SNP (effect allele)` = paste0(RSID, " (", EA, ")"),
    Cytoband = ifelse(is.na(Region) | Region == "", "NA", Region),
    `Top prioritized gene` = Top_Prioritized_Gene,
    `P-value` = fmt_p(P_meta),
    `OR (95% CI)` = fmt_or_ci(BETA_meta, SE_meta),
    `Heterogeneity P-value` = fmt_p(HetP_final),
    `I2 (%)` = ifelse(is.na(I2_final), "NA", sprintf("%.1f", I2_final)),
    coloc_text = case_when(
      !is.na(Best_PP4) & Best_PP4 >= 0.8 &
        !is.na(Coloc_Tissue) & Coloc_Tissue != "" ~
          paste0(
            "eQTL colocalization in ",
            gsub("_", " ", Coloc_Tissue),
            " (PP4=", sprintf("%.2f", Best_PP4), ")"
          ),
      TRUE ~ NA_character_
    ),
    `Biological evidence` = case_when(
      !is.na(Therapy_Evidence) & Therapy_Evidence != "" ~ Therapy_Evidence,
      !is.na(coloc_text) & coloc_text != "" ~ coloc_text,
      TRUE ~ "NA"
    )
  ) %>%
  select(
    `chr:pos:A1:A2`,
    `Lead SNP (effect allele)`,
    Cytoband,
    `Top prioritized gene`,
    Novelty,
    `P-value`,
    `OR (95% CI)`,
    `Heterogeneity P-value`,
    `I2 (%)`,
    `Biological evidence`
  )

write_tsv(out, "results/table1_compact_lead_loci.tsv", na = "NA")

md_header <- paste(names(out), collapse = " | ")
md_rule <- paste(rep("---", ncol(out)), collapse = " | ")
md_rows <- apply(out, 1, function(x) paste(x, collapse = " | "))
writeLines(
  c(
    "# Compact lead-locus summary table",
    "",
    "This version follows the short GWAS-style layout and adds one prioritized-gene column.",
    "",
    md_header,
    md_rule,
    md_rows
  ),
  "results/table1_compact_lead_loci.md"
)

cat("Wrote results/table1_compact_lead_loci.tsv\n")
cat("Wrote results/table1_compact_lead_loci.md\n")
