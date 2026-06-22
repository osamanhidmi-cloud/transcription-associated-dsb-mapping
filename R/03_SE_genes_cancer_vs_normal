# ============================================================
# 04_SE_genes_cancer_vs_normal.R
#
# Analysis module:
# SE genes in cancer vs normal
#
# Purpose:
# Compare DNA break density at super-enhancer-regulated genes
# across cancer and non-cancer cell lines.
#
# For each cell line, the break density at SE-regulated genes
# is represented relative to the median break density of non-SE
# genes from the same cell line.
#
# Input:
#   data/example/figure1d_se_gene_break_enrichment.csv
#
# Outputs:
#   results/SE_genes_cancer_vs_normal_summary.csv
#   figures/SE_genes_cancer_vs_normal_boxplot.pdf
#   figures/SE_genes_cancer_vs_normal_boxplot.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

input_file <- "data/example/figure1d_se_gene_break_enrichment.csv"
results_dir <- "results"
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Read input
# ------------------------------------------------------------

se_breaks <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "cell_line",
  "cell_type",
  "gene",
  "break_density",
  "ratio_non_se"
)

missing_cols <- setdiff(required_cols, colnames(se_breaks))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Clean and order data
# ------------------------------------------------------------

cell_line_order <- c(
  "HEK293",
  "A4098",
  "HMLE",
  "MCF10A",
  "MCF7",
  "T47D",
  "MDA-MB-436",
  "MDA-MB-468",
  "NSC",
  "SH-SY5Y"
)

se_breaks <- se_breaks |>
  mutate(
    cell_line = factor(cell_line, levels = cell_line_order),
    cell_type = factor(
      cell_type,
      levels = c("normal", "premalignant", "cancer")
    ),
    ratio_non_se = as.numeric(ratio_non_se)
  ) |>
  filter(
    !is.na(cell_line),
    !is.na(cell_type),
    !is.na(gene),
    !is.na(ratio_non_se)
  )

# ------------------------------------------------------------
# Summary per cell line
# ------------------------------------------------------------

cell_line_summary <- se_breaks |>
  group_by(cell_line, cell_type) |>
  summarise(
    n_SE_genes = n(),
    median_ratio_non_se = median(ratio_non_se, na.rm = TRUE),
    mean_ratio_non_se = mean(ratio_non_se, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cell_line_summary,
  file.path(results_dir, "SE_genes_cancer_vs_normal_summary.csv")
)

# ------------------------------------------------------------
# Statistical comparison:
# cancer versus non-cancer
# ------------------------------------------------------------

se_breaks <- se_breaks |>
  mutate(
    broad_group = ifelse(cell_type == "cancer", "cancer", "non-cancer")
  )

wilcox_cancer_vs_non_cancer <- wilcox.test(
  ratio_non_se ~ broad_group,
  data = se_breaks
)

stats_summary <- tibble(
  comparison = "Cancer versus non-cancer cell lines",
  test = "Wilcoxon rank-sum test",
  p_value = wilcox_cancer_vs_non_cancer$p.value
)

write_csv(
  stats_summary,
  file.path(results_dir, "SE_genes_cancer_vs_normal_statistics.csv")
)

# ------------------------------------------------------------
# Plot 1:
# Cell-line-level boxplot
# ------------------------------------------------------------

p_cell_lines <- ggplot(
  se_breaks,
  aes(x = cell_line, y = ratio_non_se, fill = cell_type)
) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 0.5
  ) +
  coord_cartesian(ylim = c(0, 18)) +
  labs(
    title = "DNA break enrichment at SE-regulated genes",
    subtitle = "SE genes shown relative to median non-SE genes in each cell line",
    x = "Cell line",
    y = "Break density at SE genes / median non-SE genes",
    fill = "Cell type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(figures_dir, "SE_genes_cancer_vs_normal_boxplot.pdf"),
  plot = p_cell_lines,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "SE_genes_cancer_vs_normal_boxplot.png"),
  plot = p_cell_lines,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 2:
# Broad cancer versus non-cancer comparison
# ------------------------------------------------------------

p_broad <- ggplot(
  se_breaks,
  aes(x = broad_group, y = ratio_non_se, fill = broad_group)
) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 0.5
  ) +
  coord_cartesian(ylim = c(0, 18)) +
  labs(
    title = "SE-gene break enrichment in cancer versus non-cancer cells",
    x = "",
    y = "Break density at SE genes / median non-SE genes",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(figures_dir, "SE_genes_cancer_vs_normal_grouped_boxplot.pdf"),
  plot = p_broad,
  width = 5,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "SE_genes_cancer_vs_normal_grouped_boxplot.png"),
  plot = p_broad,
  width = 5,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

message("SE genes in cancer vs normal analysis complete.")
message("Input rows: ", nrow(se_breaks))
message("Cell lines included: ", paste(unique(se_breaks$cell_line), collapse = ", "))
message("Wilcoxon cancer vs non-cancer P value: ", signif(wilcox_cancer_vs_non_cancer$p.value, 3))
message("Results written to: ", results_dir)
message("Figures written to: ", figures_dir)
