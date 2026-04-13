library(vroom)
library(dplyr)
library(data.table)

# Set working directory to project root
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# Define file paths
path_discovery <- "results/discovery_meta_results.tsv"
path_spanish <- "data/raw/Julia_2018_Spain_Only.txt"
output_path <- "results/spanish_replication_results.tsv"

# 1. Load Discovery Results
message("Loading Discovery Results...")
discovery <- vroom(path_discovery)

# 2. Extract Lead SNPs (Efficient Pruning)
message("Extracting lead SNPs (P < 5e-8)...")
# Filter by significance
sig_snps <- discovery %>% 
  filter(P_meta < 5e-8) %>%
  arrange(P_meta) %>%
  as.data.table()

if (nrow(sig_snps) == 0) {
  stop("No genome-wide significant SNPs found in discovery.")
}

message(paste("Found", nrow(sig_snps), "significant SNPs. Pruning..."))

# Fast physical distance pruning
# Use a windowed approach to avoid O(N^2) complexity
pruned_snps <- list()
while(nrow(sig_snps) > 0) {
  # Select the top SNP
  top <- sig_snps[1, ]
  pruned_snps[[length(pruned_snps) + 1]] <- top
  
  # Efficiently remove SNPs in the same region
  chr_val <- top$CHR
  bp_val <- top$BP
  window <- 500000
  
  sig_snps <- sig_snps[!(CHR == chr_val & BP > bp_val - window & BP < bp_val + window)]
}

lead_snps <- rbindlist(pruned_snps)
message(paste("Found", nrow(lead_snps), "independent lead loci."))

# 3. Load Spanish Replication Data
message("Loading Spanish Replication Data...")
# Columns: CHR, SNP (RSID), BP, A1leleA, AlleleB, OR, OR_lower, OR_upper, P
# Verification shows AlleleB is likely the Effect Allele (A1) in this dataset.
spanish <- vroom(path_spanish) %>%
  select(RSID = SNP, OA_rep = A1leleA, EA_rep = AlleleB, OR_rep = OR, P_rep = P) %>%
  mutate(BETA_rep = log(OR_rep))

# 4. Join and Align
message("Aligning with replication data...")
replicated_loci <- inner_join(lead_snps, spanish, by = "RSID")

# Alignment function
align_to_meta <- function(df) {
  # Align BETA_rep to EA (Effect Allele from Meta)
  # OA OA_rep, EA EA_rep
  
  idx_same <- df$EA == df$EA_rep & df$OA == df$OA_rep
  idx_flip <- df$EA == df$OA_rep & df$OA == df$EA_rep
  
  # Filter to matched
  df <- df[idx_same | idx_flip, ]
  idx_flip_new <- df$EA == df$OA_rep
  
  df$BETA_rep[idx_flip_new] <- -df$BETA_rep[idx_flip_new]
  
  # Concordance check
  df$concordant <- sign(df$BETA_meta) == sign(df$BETA_rep)
  df$replicated <- df$P_rep < 0.05 & df$concordant
  
  return(df)
}

results <- align_to_meta(replicated_loci)
message(paste("Replicated", sum(results$replicated), "of", nrow(results), "loci."))

# 5. Save Results
vroom_write(results, output_path)
message("Done.")
