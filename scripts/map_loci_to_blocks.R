# scripts/map_loci_to_blocks.R
library(data.table)

# Paths
TOP_LOCI_FILE <- "results/top_loci_summary_table.tsv"
BLOCKS_FILE <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/uterine_fibroids/mr_analysis_2026_04_16/reference_data/LAVA_s2500_m25_f1_w200.blocks"

# Load
top <- fread(TOP_LOCI_FILE)
blocks <- fread(BLOCKS_FILE)
blocks[, LOC := 1:.N]

# Map
top[, LAVA_LOC := as.integer(NA)]

for (i in 1:nrow(top)) {
  CHR_val <- top$CHR[i]
  BP_val <- top$BP[i]
  
  match <- blocks[chr == CHR_val & start <= BP_val & stop >= BP_val]
  if (nrow(match) > 0) {
    top$LAVA_LOC[i] <- match$LOC[1]
  }
}

# Save mapping
fwrite(top, "results/top_loci_lava_mapping.tsv", sep="\t")
cat("Mapped", sum(!is.na(top$LAVA_LOC)), "out of", nrow(top), "loci to LAVA blocks.\n")
Unique_LOCs <- unique(top$LAVA_LOC[!is.na(top$LAVA_LOC)])
cat("These loci cover", length(Unique_LOCs), "unique LAVA blocks.\n")
cat("Critical Blocks to watch:", paste(sort(Unique_LOCs), collapse=", "), "\n")
