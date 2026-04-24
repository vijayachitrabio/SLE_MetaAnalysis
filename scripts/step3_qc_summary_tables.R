library(vroom)
library(dplyr)
library(data.table)


# Set working directory to project root

# Define file paths
path_rep <- "results/spanish_replication_results.tsv"
path_anno <- "results/gene_annotation.tsv"
output_table <- "results/top_loci_summary_table.tsv"
output_qc <- "results/supplementary_qc_table.tsv"

# 1. Load Data
message("Loading replication results and annotations...")
rep_results <- vroom(path_rep)

# Load existing annotations (RSID -> Gene/Region mapping)
if (file.exists(path_anno)) {
  anno <- vroom(path_anno) %>%
    select(RSID = SNP, Gene, Region, Known_SLE) %>%
    distinct()
} else {
  anno <- data.frame(RSID = character(), Gene = character(), Region = character(), Known_SLE = logical())
}

# 2. Merge and Format
message("Formatting final table...")
final_table <- rep_results %>%
  left_join(anno, by = "RSID") %>%
  mutate(
    Known_SLE = ifelse(is.na(Known_SLE), FALSE, Known_SLE),
    Gene = ifelse(is.na(Gene), "TBD", Gene),
    Region = ifelse(is.na(Region), "TBD", Region)
  ) %>%
  select(
    RSID, CHR, BP, OA, EA, 
    BETA_disco = BETA_meta, P_disco = P_meta,
    BETA_rep, P_rep, 
    I2, HetP,
    Gene, Region, Known_SLE, Replicated = replicated
  )

# 3. Handle HLA special note
final_table <- final_table %>%
  mutate(
    Notes = ifelse(CHR == 6 & BP > 25000000 & BP < 35000000, 
                   "HLA region; High I2 expected due to North-South LD divergence", "")
  )

# 4. Save Final Tables
message("Saving final tables...")
vroom_write(final_table, output_table)

# Supplementary QC Table (Discovery Details)
qc_table <- final_table %>%
  left_join(rep_results %>% select(RSID, P_bentham, P_finngen), by = "RSID") %>%
  select(RSID, CHR, BP, P_disco, P_bentham, P_finngen, I2, HetP, Notes)

vroom_write(qc_table, output_qc)
message("Done.")
