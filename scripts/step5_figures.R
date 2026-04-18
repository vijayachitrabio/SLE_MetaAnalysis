library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")
dir.create("figures", showWarnings = FALSE)

message("Loading full discovery meta-results for Manhattan/QQ...")
# Load full data
meta_data <- fread("results/discovery_meta_results.tsv", select = c("CHR", "BP", "P_meta", "RSID"))

# Only keep snps with reasonable p-values for plotting to save memory
plot_data <- meta_data[P_meta < 0.05]
rm(meta_data)
gc()

# Calculate cumulative BP for Manhattan
plot_data <- plot_data %>% 
  arrange(CHR, BP) %>%
  mutate(CHR = as.numeric(CHR)) %>%
  filter(!is.na(CHR))

data_cum <- plot_data %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

plot_data <- plot_data %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add) %>%
  mutate(minuslog10p = -log10(P_meta))

axisdf <- plot_data %>% group_by(CHR) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )

message("Generating Manhattan Plot...")

# Identify Top Genes for labels (High-Confidence & Replicated)
top_loci <- fread("results/top_loci_summary_table.tsv")
label_data <- top_loci %>%
  filter(Replicated == TRUE) %>%
  mutate(label = Gene) %>%
  arrange(P_meta) %>%
  head(15)

# Merge with bp_cum 
label_points <- label_data %>%
  inner_join(plot_data %>% select(RSID, bp_cum), by = "RSID")

p_manhattan <- ggplot(plot_data, aes(x=bp_cum, y=minuslog10p)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#2c3e50", "#2980b9"), 22)) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = -log10(5e-8), color = "#e74c3c", linetype="dashed") +
  # Add Labels
  ggrepel::geom_text_repel(
    data = label_points,
    aes(x = bp_cum, y = -log10(P_meta), label = label),
    size = 3.2, fontface = "italic", color = "black",
    box.padding = 0.5, point.padding = 0.5, force = 4,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(x = "Chromosome", y = "-log10(P-value)", title = "North-to-South European SLE Meta-Analysis")

ggsave("figures/manhattan_plot.png", p_manhattan, width = 14, height = 7, dpi = 300)

message("Generating QQ Plot...")
# Downsample for QQ plot speed
qq_data <- plot_data %>% sample_frac(0.1)
y <- -log10(sort(qq_data$P_meta, decreasing = FALSE))
x <- -log10(ppoints(length(y)))

p_qq <- ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) +
  geom_point(color="#2980b9", alpha=0.5, size=1) +
  geom_abline(intercept = 0, slope = 1, color="#e74c3c") +
  theme_minimal() +
  labs(x="Expected -log10(P)", y="Observed -log10(P)", title="QQ Plot of Meta-Analysis")

ggsave("figures/qq_plot.png", p_qq, width = 6, height = 6, dpi = 300)

message("Generating Forest Plot for Replicated Loci...")
# Load top loci and replication details
top_loci <- fread("results/top_loci_summary_table.tsv")
rep_results <- fread("results/spanish_replication_results.tsv")

# Filter for High-Confidence Replicated Loci
rep_loci <- top_loci[Replicated == TRUE][order(P_meta)][head(1:nrow(top_loci), 15)]

# Merge with detailed replication stats to get BETA_rep
# Note: P_rep is already in top_loci from step21, so we exclude it from the merge to avoid .x/.y suffixes
forest_data <- merge(rep_loci, rep_results[, .(RSID, BETA_rep, SE_rep_manual = SE_meta)], by = "RSID", all.x = TRUE)

# Calculate Discovery SE from BETA and P_meta
forest_data[, SE_disco := abs(BETA / qnorm(P_meta/2))]
# Calculate Replication SE from BETA_rep and P_rep
forest_data[, SE_rep := abs(BETA_rep / qnorm(P_rep/2))]

# Prepare for plotting
forest_df <- bind_rows(
  forest_data %>% select(RSID, Gene, BETA_val = BETA, SE_val = SE_disco) %>% mutate(Cohort = "Discovery (North)"),
  forest_data %>% select(RSID, Gene, BETA_val = BETA_rep, SE_val = SE_rep) %>% mutate(Cohort = "Replication (South)")
) %>%
  mutate(Label = ifelse(Gene == "" | is.na(Gene), RSID, paste(Gene, "\n", RSID)))

p_forest <- ggplot(forest_df, aes(x=BETA_val, y=Label, color=Cohort)) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbarh(aes(xmin=BETA_val - 1.96*SE_val, xmax=BETA_val + 1.96*SE_val), position=position_dodge(width=0.5), height=0.2) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  theme_minimal() +
  labs(x="Effect Size (Beta)", y="", title="Effect Size Comparison (Top Replicated Loci)") +
  scale_color_manual(values=c("#2980b9", "#e67e22"))

ggsave("figures/forest_plot_top15.png", p_forest, width = 8, height = 8, dpi = 300)

message("Figures generated successfully.")
