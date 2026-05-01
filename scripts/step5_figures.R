library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd(".")


dir.create("figures", showWarnings = FALSE)

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub('^"|"$', "", x))
  x[x %in% c("", "NA", "TBD", "character(0)", "Locus_NA")] <- NA_character_
  x
}

extract_near_gene <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- grepl("Known \\(near .+\\)", x)
  out[hit] <- sub("^Known \\(near (.+)\\)$", "\\1", x[hit])
  out
}

choose_gene_label <- function(gene, master_gene, effector_gene, novelty, rsid) {
  label <- dplyr::coalesce(
    clean_text(gene),
    clean_text(master_gene),
    clean_text(effector_gene),
    extract_near_gene(as.character(novelty))
  )
  label <- ifelse(grepl(" / ", label), sub("^.* / ", "", label), label)
  ifelse(is.na(label), rsid, label)
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

message("Loading full discovery meta-results for Manhattan/QQ...")
# Load full data
meta_data <- fread("results/discovery_meta_results.tsv", select = c("CHR", "BP", "P_meta", "RSID"))

# Calculate true Lambda GC on the FULL dataset before filtering
message("Calculating genome-wide Lambda GC...")
lambda_gc <- median(qchisq(meta_data$P_meta, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
message("Calculated Meta-Analysis Lambda: ", round(lambda_gc, 3))

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
master_loci <- fread("results/master_results_table.tsv")

master_labels <- master_loci %>%
  mutate(
    Master_Gene = clean_text(Gene),
    Master_Effector_Gene = clean_text(Effector_Gene)
  ) %>%
  group_by(RSID) %>%
  summarise(
    Master_Gene = first_non_missing(Master_Gene),
    Master_Effector_Gene = first_non_missing(Master_Effector_Gene),
    Novelty = dplyr::first(Novelty),
    .groups = "drop"
  )

label_data <- top_loci %>%
  left_join(master_labels, by = "RSID") %>%
  filter(Replicated == TRUE) %>%
  mutate(label = choose_gene_label(Gene, Master_Gene, Master_Effector_Gene, Novelty.x, RSID)) %>%
  arrange(P_meta) %>%
  head(7)

# Merge with bp_cum 
label_points <- label_data %>%
  inner_join(plot_data %>% select(RSID, bp_cum) %>% distinct(), by = "RSID")

clic1_label <- label_points %>% filter(label == "CLIC1")
other_labels <- label_points %>% filter(label != "CLIC1")

p_manhattan <- ggplot(plot_data, aes(x=bp_cum, y=minuslog10p)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#2c3e50", "#2980b9"), 22)) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  geom_hline(yintercept = -log10(5e-8), color = "#e74c3c", linetype="dashed") +
  # Add Labels
  ggrepel::geom_text_repel(
    data = other_labels,
    aes(x = bp_cum, y = -log10(P_meta), label = label),
    size = 3.2, fontface = "plain", color = "black",
    box.padding = 0.5, point.padding = 0.5, force = 4,
    max.overlaps = Inf, min.segment.length = 0, seed = 123
  ) +
  ggrepel::geom_text_repel(
    data = clic1_label,
    aes(x = bp_cum, y = -log10(P_meta), label = label),
    size = 3.5, fontface = "bold", color = "black",
    box.padding = 0.6, point.padding = 0.6, force = 5,
    max.overlaps = Inf, min.segment.length = 0, seed = 123
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.margin = margin(12, 20, 12, 12)) +
  labs(x = "Chromosome", y = "-log10(P-value)", title = NULL)

ggsave("figures/manhattan_plot.png", p_manhattan, width = 14, height = 7, dpi = 300)
ggsave("figures/manhattan_plot.pdf", p_manhattan, width = 14, height = 7)

message("Generating QQ Plots (Before and After GC)...")

# Downsample for QQ plot speed
qq_data <- plot_data %>% sample_frac(0.1)

# 1. Before GC (Raw P_meta)
y_before <- -log10(sort(qq_data$P_meta, decreasing = FALSE))
x_before <- -log10(ppoints(length(y_before)))
p_qq_before <- ggplot(data.frame(x=x_before, y=y_before), aes(x=x, y=y)) +
  geom_point(color="#2980b9", alpha=0.5, size=1) +
  geom_abline(intercept = 0, slope = 1, color="#e74c3c") +
  theme_minimal() +
  labs(x="Expected -log10(P)", y="Observed -log10(P)", 
       title = sprintf("Before GC (\u03bb = %.2f)", lambda_gc))

# 2. After GC (Corrected P_meta)
# Apply Genomic Control by scaling the chi-squared statistic
qq_data <- qq_data %>%
  mutate(chisq_raw = qchisq(P_meta, df = 1, lower.tail = FALSE),
         chisq_gc = chisq_raw / lambda_gc,
         P_gc = pchisq(chisq_gc, df = 1, lower.tail = FALSE))

y_after <- -log10(sort(qq_data$P_gc, decreasing = FALSE))
x_after <- -log10(ppoints(length(y_after)))
p_qq_after <- ggplot(data.frame(x=x_after, y=y_after), aes(x=x, y=y)) +
  geom_point(color="#27ae60", alpha=0.5, size=1) +
  geom_abline(intercept = 0, slope = 1, color="#e74c3c") +
  theme_minimal() +
  labs(x="Expected -log10(P)", y="Observed -log10(P)", 
       title = "After GC (\u03bb = 1.00)")

# Combine plots side-by-side
p_qq_combined <- cowplot::plot_grid(p_qq_before, p_qq_after, labels = c("A", "B"))
ggsave("figures/qq_plot_combined.png", p_qq_combined, width = 10, height = 5, dpi = 300)
message("Saved Before/After GC QQ plots to figures/qq_plot_combined.png")

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
