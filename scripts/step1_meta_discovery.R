library(vroom)
library(dplyr)
library(data.table)

# Set working directory to project root
setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")

# Define file paths
path_bentham <- "data/raw/Bentham_2015_SLE.h.tsv.gz"
path_finngen <- "data/raw/finngen_R12_M13_SLE.gz"
output_path <- "results/discovery_meta_results.tsv"

# 1. Load Bentham (Discovery 1 - Northern European)
message("Loading Bentham (Discovery 1)...")
# Columns: chromosome, base_pair_location, other_allele, effect_allele, beta, standard_error, p_value, hm_rsid
bentham <- vroom(path_bentham, col_select = c(chromosome, base_pair_location, other_allele, effect_allele, beta, standard_error, p_value, hm_rsid)) %>%
  rename(CHR = chromosome, BP = base_pair_location, OA = other_allele, EA = effect_allele, BETA = beta, SE = standard_error, P = p_value) %>%
  filter(!is.na(BETA), !is.na(SE)) %>%
  mutate(variant_id = paste(CHR, BP, sep = ":"))

# 2. Load FinnGen (Discovery 2 - North-Eastern European)
message("Loading FinnGen (Discovery 2)...")
# Columns: #chrom, pos, ref, alt, beta, sebeta, pval
finngen <- vroom(path_finngen, col_select = c(`#chrom`, pos, ref, alt, beta, sebeta, pval)) %>%
  rename(CHR = `#chrom`, BP = pos, OA = ref, EA = alt, BETA = beta, SE = sebeta, P = pval) %>%
  filter(!is.na(BETA), !is.na(SE)) %>%
  mutate(variant_id = paste(CHR, BP, sep = ":"))

# 3. Harmonize Alleles (Align FinnGen to Bentham)
message("Aligning alleles...")
merged <- inner_join(bentham, finngen, by = "variant_id", suffix = c("_b", "_f"))

# Function to align effect alleles
align_alleles <- function(df) {
  # df has BETA_b, SE_b, OA_b, EA_b and BETA_f, SE_f, OA_f, EA_f
  # We align everything to EA_b (Bentham's Effect Allele)
  
  # Same alleles
  idx_same <- df$EA_b == df$EA_f & df$OA_b == df$OA_f
  # Flipped alleles
  idx_flip <- df$EA_b == df$OA_f & df$OA_b == df$EA_f
  
  # Filter to matched variants
  df <- df[idx_same | idx_flip, ]
  idx_flip_new <- df$EA_b == df$OA_f
  
  # Flip FinnGen effects where necessary
  df$BETA_f[idx_flip_new] <- -df$BETA_f[idx_flip_new]
  
  return(df)
}

harmonized <- align_alleles(merged)
message(paste("Harmonized", nrow(harmonized), "variants."))

# 4. Meta-Analysis (Inverse Variance Weighted Fixed Effects)
message("Running inverse-variance meta-analysis...")
w_b <- 1 / (harmonized$SE_b^2)
w_f <- 1 / (harmonized$SE_f^2)

# Meta Beta
harmonized$meta_beta <- (harmonized$BETA_b * w_b + harmonized$BETA_f * w_f) / (w_b + w_f)
# Meta SE
harmonized$meta_se <- sqrt(1 / (w_b + w_f))
# Meta P-value (Z-test)
harmonized$meta_z <- harmonized$meta_beta / harmonized$meta_se
harmonized$meta_p <- 2 * pnorm(-abs(harmonized$meta_z))

# 5. Heterogeneity Stats
# Cochran's Q
harmonized$Q <- w_b * (harmonized$BETA_b - harmonized$meta_beta)^2 + w_f * (harmonized$BETA_f - harmonized$meta_beta)^2
# df = k - 1 = 1
harmonized$HetP <- pchisq(harmonized$Q, df = 1, lower.tail = FALSE)
# I^2 = max(0, (Q - df) / Q)
harmonized$I2 <- pmax(0, (harmonized$Q - 1) / harmonized$Q) * 100

# 6. Save Results
message("Saving results...")
final_output <- harmonized %>%
  select(CHR = CHR_b, BP = BP_b, RSID = hm_rsid, OA = OA_b, EA = EA_b, 
         BETA_bentham = BETA_b, SE_bentham = SE_b, P_bentham = P_b,
         BETA_finngen = BETA_f, SE_finngen = SE_f, P_finngen = P_f,
         BETA_meta = meta_beta, SE_meta = meta_se, P_meta = meta_p,
         I2, HetP)

vroom_write(final_output, output_path)
message("Done.")
