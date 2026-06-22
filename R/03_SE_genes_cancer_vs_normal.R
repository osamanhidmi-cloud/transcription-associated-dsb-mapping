# ============================================================
# 03_SE_genes_cancer_vs_normal.R
#
# Analysis module:
# SE genes in cancer vs normal
#
# Purpose:
# Generate the Figure 1D-style boxplot:
# break density at SE-regulated genes relative to all non-SE genes
# across multiple cell lines.
#
# This script intentionally uses base R boxplot to stay close to
# the original plotting code used for the paper.
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
# Clean data
# ------------------------------------------------------------

se_breaks <- se_breaks |>
  mutate(
    cell_line = as.character(cell_line),
    cell_type = as.character(cell_type),
    gene = as.character(gene),
    ratio_non_se = as.numeric(ratio_non_se)
  ) |>
  filter(
    !is.na(cell_line),
    !is.na(gene),
    !is.na(ratio_non_se)
  )

# ------------------------------------------------------------
# Cell-line order exactly matching the original boxplot order
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

# ------------------------------------------------------------
# Create vectors for each cell line
# ------------------------------------------------------------

HEK293_ratio <- se_breaks |>
  filter(cell_line == "HEK293") |>
  pull(ratio_non_se)

A4098_ratio <- se_breaks |>
  filter(cell_line == "A4098") |>
  pull(ratio_non_se)

HMLE_ratio <- se_breaks |>
  filter(cell_line == "HMLE") |>
  pull(ratio_non_se)

MCF10A_ratio <- se_breaks |>
  filter(cell_line == "MCF10A") |>
  pull(ratio_non_se)

MCF7_ratio <- se_breaks |>
  filter(cell_line == "MCF7") |>
  pull(ratio_non_se)

T47D_ratio <- se_breaks |>
  filter(cell_line == "T47D") |>
  pull(ratio_non_se)

MDA_MB_436_ratio <- se_breaks |>
  filter(cell_line == "MDA-MB-436") |>
  pull(ratio_non_se)

MDA_MB_468_ratio <- se_breaks |>
  filter(cell_line == "MDA-MB-468") |>
  pull(ratio_non_se)

NSC_ratio <- se_breaks |>
  filter(cell_line == "NSC") |>
  pull(ratio_non_se)

SH_SY5Y_ratio <- se_breaks |>
  filter(cell_line == "SH-SY5Y") |>
  pull(ratio_non_se)

# ------------------------------------------------------------
# Export simple summary table
# ------------------------------------------------------------

summary_tbl <- se_breaks |>
  mutate(cell_line = factor(cell_line, levels = cell_line_order)) |>
  group_by(cell_line, cell_type) |>
  summarise(
    n_SE_genes = n(),
    median_ratio_non_se = median(ratio_non_se, na.rm = TRUE),
    mean_ratio_non_se = mean(ratio_non_se, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(cell_line)

write_csv(
  summary_tbl,
  file.path(results_dir, "SE_genes_cancer_vs_normal_summary.csv")
)

# ------------------------------------------------------------
# Base R plotting function
# Close to the original Figure 1D plotting code
# ------------------------------------------------------------

plot_se_boxplot <- function() {
  
  boxplot(
    HEK293_ratio,
    A4098_ratio,
    HMLE_ratio,
    MCF10A_ratio,
    MCF7_ratio,
    T47D_ratio,
    MDA_MB_436_ratio,
    MDA_MB_468_ratio,
    NSC_ratio,
    SH_SY5Y_ratio,
    ylab = "Break density (SE genes relative to all non-SE genes)",
    cex.lab = 1.4,
    outpch = 16,
    outcex = 0.2,
    col = c(
      "#dd9999",
      "#aa1133",
      "gray",
      "#dd9999",
      "#aa1133",
      "#aa1133",
      "#aa1133",
      "#aa1133",
      "gray",
      "#aa1133"
    ),
    ylim = c(0, 18)
  )
  
  y_ticks <- c(seq(0, 5, by = 1), seq(10, 15, by = 5))
  y_labels <- y_ticks
  
  axis(2, at = y_ticks, labels = y_labels)
}

# ------------------------------------------------------------
# Save PDF
# ------------------------------------------------------------

pdf(
  file = file.path(figures_dir, "SE_genes_cancer_vs_normal_boxplot.pdf"),
  width = 8,
  height = 5
)

plot_se_boxplot()

dev.off()

# ------------------------------------------------------------
# Save PNG
# ------------------------------------------------------------

png(
  filename = file.path(figures_dir, "SE_genes_cancer_vs_normal_boxplot.png"),
  width = 8,
  height = 5,
  units = "in",
  res = 300
)

plot_se_boxplot()

dev.off()

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

message("SE genes in cancer vs normal boxplot complete.")
message("Input rows: ", nrow(se_breaks))
message("Summary written to: ", file.path(results_dir, "SE_genes_cancer_vs_normal_summary.csv"))
message("Figure written to: ", figures_dir)
