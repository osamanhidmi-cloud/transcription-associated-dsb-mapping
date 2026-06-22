# ============================================================
# 03_plot_transcription_stress_results.R
#
# Purpose:
# Plot transcription-stress scores and gene-set enrichment results.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

results_dir <- "results"
figures_dir <- "figures"

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

n_genes_to_plot <- 30

# ------------------------------------------------------------
# Read results
# ------------------------------------------------------------

tss_scores <- read_csv(
  file.path(results_dir, "transcription_stress_gene_scores.csv"),
  show_col_types = FALSE
)

enrichment_results <- read_csv(
  file.path(results_dir, "transcription_stress_gene_set_enrichment.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Plot 1: ranked transcription-stress genes
# ------------------------------------------------------------

top_tss_plot <- tss_scores |>
  arrange(desc(TSS_final)) |>
  slice_head(n = n_genes_to_plot) |>
  mutate(
    gene = factor(gene, levels = gene),
    SE_regulated = ifelse(gene %in% gene[SE_regulated %in% TRUE], TRUE, FALSE)
  )

p_tss <- ggplot(
  top_tss_plot,
  aes(x = gene, y = TSS_final, fill = SE_regulated)
) +
  geom_col() +
  labs(
    title = "Top transcription-stress genes",
    x = "Gene",
    y = "Final transcription-stress score",
    fill = "SE-regulated"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#cc0000", "FALSE" = "#0000cc")
  )

ggsave(
  filename = file.path(figures_dir, "top_transcription_stress_genes_barplot.pdf"),
  plot = p_tss,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "top_transcription_stress_genes_barplot.png"),
  plot = p_tss,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 2: heatmap of individual z-scores
# ------------------------------------------------------------

heatmap_df <- tss_scores |>
  arrange(desc(TSS_final)) |>
  slice_head(n = n_genes_to_plot) |>
  select(gene, dsb_z, rloop_z, top1cc_z, top1_z)

heatmap_mat <- heatmap_df |>
  tibble::column_to_rownames("gene") |>
  as.matrix()

pdf(file.path(figures_dir, "top_transcription_stress_genes_heatmap.pdf"),
    width = 5,
    height = 8)

pheatmap(
  heatmap_mat,
  color = colorRampPalette(c("#0000cc", "white", "#cc0000"))(100),
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "white",
  fontsize_row = 7,
  fontsize_col = 9
)

dev.off()

# ------------------------------------------------------------
# Plot 3: enrichment ratio plot
# ------------------------------------------------------------

enrichment_plot_df <- enrichment_results |>
  filter(direction == "enrichment") |>
  arrange(desc(enrichment_ratio)) |>
  mutate(
    gene_set = factor(gene_set, levels = gene_set),
    minus_log10_p_capped = pmin(minus_log10_p, 30)
  )

p_enrichment <- ggplot(
  enrichment_plot_df,
  aes(
    x = gene_set,
    y = enrichment_ratio + 0.1,
    fill = minus_log10_p_capped
  )
) +
  geom_col() +
  scale_y_log10() +
  labs(
    title = "Gene-set enrichment among transcription-stress genes",
    x = "Gene set",
    y = "Observed / expected overlap",
    fill = "-log10(P value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_gradient2(
    low = "#0000cc",
    mid = "#cc0000",
    high = "#cc0000",
    midpoint = 20,
    limits = c(0, 30)
  )

ggsave(
  filename = file.path(figures_dir, "transcription_stress_gene_set_enrichment.pdf"),
  plot = p_enrichment,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "transcription_stress_gene_set_enrichment.png"),
  plot = p_enrichment,
  width = 8,
  height = 5,
  dpi = 300
)

message("Plots written to: ", figures_dir)
