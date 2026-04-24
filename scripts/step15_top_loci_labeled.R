#!/usr/bin/env Rscript
# Top loci figure with improved gene labels and ranked visualization

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})


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
  ifelse(is.na(label), rsid, label)
}

loci <- fread("results/top_loci_summary_table.tsv")
master <- fread("results/master_results_table.tsv")

master_best <- master %>%
  mutate(
    Gene = clean_text(Gene),
    Effector_Gene = clean_text(Effector_Gene),
    has_effector = !is.na(Effector_Gene),
    has_gene = !is.na(Gene)
  ) %>%
  arrange(desc(has_effector), desc(has_gene)) %>%
  group_by(RSID) %>%
  slice(1) %>%
  ungroup() %>%
  select(RSID, Master_Gene = Gene, Master_Effector_Gene = Effector_Gene)

top_loci <- loci %>%
  left_join(master_best, by = "RSID") %>%
  mutate(
    Plot_Gene = choose_gene_label(Gene, Master_Gene, Master_Effector_Gene, Novelty, RSID),
    Plot_Gene = ifelse(grepl(" / ", Plot_Gene), sub("^.* / ", "", Plot_Gene), Plot_Gene),
    logP = -log10(P_meta),
    abs_beta = abs(BETA),
    Direction = ifelse(BETA >= 0, "Risk increasing", "Protective"),
    Replication_Status = ifelse(Replicated, "LAVA replicated", "Not replicated"),
    Novelty_Group = ifelse(grepl("^Putative Novel$", Novelty), "Putative novel", "Known / nearby known"),
    Label = paste0(Plot_Gene, "  ", RSID)
  ) %>%
  arrange(P_meta, desc(abs_beta)) %>%
  slice_head(n = 10) %>%
  mutate(Label = factor(Label, levels = rev(Label)))

subtitle_text <- paste(
  "Top 10 loci ranked by discovery meta-analysis significance.",
  "Color shows novelty, shape shows LAVA replication, point size shows |beta|."
)

p_dot <- ggplot(top_loci, aes(x = logP, y = Label)) +
  geom_segment(aes(x = 0, xend = logP, y = Label, yend = Label), linewidth = 0.8, color = "#d7dde5") +
  geom_point(
    aes(size = abs_beta, fill = Novelty_Group, shape = Replication_Status),
    color = "#1f2933", stroke = 0.5
  ) +
  scale_fill_manual(values = c(
    "Putative novel" = "#d95f02",
    "Known / nearby known" = "#1b9e77"
  )) +
  scale_shape_manual(values = c(
    "LAVA replicated" = 21,
    "Not replicated" = 24
  )) +
  scale_size_continuous(name = "|Beta|", range = c(4, 9)) +
  labs(
    title = "Top 10 SLE-Associated Loci",
    subtitle = subtitle_text,
    x = expression(-log[10](P[meta])),
    y = NULL,
    fill = "Locus class",
    shape = "Support"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "#18212b"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, color = "#4a5565")
  )

p_lollipop <- ggplot(top_loci, aes(x = logP, y = Label)) +
  geom_segment(aes(x = 0, xend = logP, y = Label, yend = Label, color = Novelty_Group), linewidth = 1.1, alpha = 0.75) +
  geom_point(
    aes(size = abs_beta, fill = Novelty_Group, shape = Replication_Status),
    color = "#18212b", stroke = 0.5
  ) +
  geom_text(
    aes(
      x = pmin(logP + max(logP) * 0.02, max(logP) * 1.08),
      label = paste0("beta=", sprintf("%.2f", BETA))
    ),
    hjust = 0, size = 3.3, color = "#46505c"
  ) +
  scale_color_manual(values = c(
    "Putative novel" = "#d95f02",
    "Known / nearby known" = "#1b9e77"
  ), guide = "none") +
  scale_fill_manual(values = c(
    "Putative novel" = "#d95f02",
    "Known / nearby known" = "#1b9e77"
  )) +
  scale_shape_manual(values = c(
    "LAVA replicated" = 21,
    "Not replicated" = 24
  )) +
  scale_size_continuous(name = "|Beta|", range = c(4.2, 8.5)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Top 10 SLE Loci Ranked by Statistical Strength",
    subtitle = "Top 10 loci ranked by discovery meta-analysis significance.",
    x = expression(-log[10](P[meta])),
    y = NULL,
    fill = "Locus class",
    shape = "Support"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "#18212b"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, color = "#4a5565")
  )

ggsave("figures/top_loci_labeled.png", p_dot, width = 11, height = 7.5, dpi = 320, bg = "white")
ggsave("figures/top_loci_labeled.pdf", p_dot, width = 11, height = 7.5)
ggsave("figures/top_loci_ranked_lollipop.png", p_lollipop, width = 11.5, height = 7.5, dpi = 320, bg = "white")
ggsave("figures/top_loci_ranked_lollipop.pdf", p_lollipop, width = 11.5, height = 7.5)

cat("Saved improved top-loci dot plot and alternative lollipop plot to figures/.\n")
print(top_loci %>% select(RSID, Plot_Gene, Novelty, Replicated, logP, BETA))
