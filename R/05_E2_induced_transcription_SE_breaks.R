# ============================================================
# 05_E2_induced_transcription_SE_breaks.R
#
# Analysis module:
# E2-induced transcription and SE-associated DSBs
#
# Purpose:
# Generate a Figure 3D-style boxplot comparing DSB density
# before and after E2 treatment across GRO-seq induction groups.
#
# Genes are grouped by E2-induced GRO-seq fold-change and split
# into:
#   1. non-ERSE / non-SE-associated genes
#   2. ERSE / SE-associated genes
#
# Input:
#   data/example/figure3d_E2_GROseq_SE_breaks.csv
#
# Outputs:
#   results/E2_induced_transcription_SE_breaks_summary.csv
#   figures/E2_induced_transcription_SE_breaks_boxplot.pdf
#   figures/E2_induced_transcription_SE_breaks_boxplot.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

input_file <- "data/example/figure3d_E2_GROseq_SE_breaks.csv"
results_dir <- "results"
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Read input
# ------------------------------------------------------------

fig3d <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "gene",
  "FC",
  "FC_group",
  "SE_status",
  "NE",
  "EE"
)

missing_cols <- setdiff(required_cols, colnames(fig3d))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Clean and order data
# ------------------------------------------------------------

fc_group_order <- c("1-2", "2-3", "3-4", "4-7", "7-10", ">10")

fig3d <- fig3d |>
  mutate(
    gene = as.character(gene),
    FC_group = factor(FC_group, levels = fc_group_order),
    SE_status = factor(SE_status, levels = c("non_ERSE", "ERSE")),
    NE = as.numeric(NE),
    EE = as.numeric(EE)
  ) |>
  filter(
    !is.na(gene),
    !is.na(FC_group),
    !is.na(SE_status),
    !is.na(NE),
    !is.na(EE)
  )

# ------------------------------------------------------------
# Export summary table
# ------------------------------------------------------------

summary_tbl <- fig3d |>
  group_by(FC_group, SE_status) |>
  summarise(
    n_genes = n(),
    median_NE = median(NE, na.rm = TRUE),
    median_EE = median(EE, na.rm = TRUE),
    mean_NE = mean(NE, na.rm = TRUE),
    mean_EE = mean(EE, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  summary_tbl,
  file.path(results_dir, "E2_induced_transcription_SE_breaks_summary.csv")
)

# ------------------------------------------------------------
# Convert table into boxplot input list
# Order:
#   non-ERSE NE, non-ERSE EE, ERSE NE, ERSE EE, gap
# repeated for each FC group
# ------------------------------------------------------------

boxplot_values <- list()
boxplot_colors <- c()
group_labels <- c()

for (grp in fc_group_order) {
  
  non_erse <- fig3d |>
    filter(FC_group == grp, SE_status == "non_ERSE")
  
  erse <- fig3d |>
    filter(FC_group == grp, SE_status == "ERSE")
  
  boxplot_values <- c(
    boxplot_values,
    list(
      non_erse$NE,
      non_erse$EE,
      erse$NE,
      erse$EE,
      NA
    )
  )
  
  boxplot_colors <- c(
    boxplot_colors,
    "#cccccc",
    "#cccccc",
    "#ee7777",
    "#ee7777",
    NA
  )
  
  group_labels <- c(
    group_labels,
    "-",
    "+",
    "-",
    "+",
    ""
  )
}

# Remove final gap
boxplot_values <- boxplot_values[-length(boxplot_values)]
boxplot_colors <- boxplot_colors[-length(boxplot_colors)]
group_labels <- group_labels[-length(group_labels)]

# ------------------------------------------------------------
# Base R plotting function
# Close to the original Figure 3D plotting code
# ------------------------------------------------------------

plot_fig3d_boxplot <- function() {
  
  bp <- boxplot(
    boxplot_values,
    col = boxplot_colors,
    ylim = c(0, 600),
    outpch = 16,
    outcex = 0.2,
    ylab = "Break density",
    whisklty = 1,
    outline = FALSE,
    xaxt = "n",
    las = 2
  )
  
  label_positions <- which(group_labels != "")
  
  axis(
    side = 1,
    at = label_positions,
    labels = group_labels[label_positions]
  )
  
  # Optional group labels below the -/+ labels
  group_centers <- seq(2.5, by = 5, length.out = length(fc_group_order))
  
  mtext(
    text = fc_group_order,
    side = 1,
    at = group_centers,
    line = 2.2,
    cex = 0.8
  )
  
  mtext(
    text = "GRO-seq E2 induction group",
    side = 1,
    line = 3.5,
    cex = 0.9
  )
  
  legend(
    "topright",
    legend = c("non-ERSE", "ERSE"),
    fill = c("#cccccc", "#ee7777"),
    border = "black",
    bty = "n",
    cex = 0.8
  )
}

# ------------------------------------------------------------
# Save PDF
# ------------------------------------------------------------

pdf(
  file = file.path(figures_dir, "E2_induced_transcription_SE_breaks_boxplot.pdf"),
  width = 9,
  height = 5
)

plot_fig3d_boxplot()

dev.off()

# ------------------------------------------------------------
# Save PNG
# ------------------------------------------------------------

png(
  filename = file.path(figures_dir, "E2_induced_transcription_SE_breaks_boxplot.png"),
  width = 9,
  height = 5,
  units = "in",
  res = 300
)

plot_fig3d_boxplot()

dev.off()

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

message("E2-induced transcription and SE-associated DSB analysis complete.")
message("Input rows: ", nrow(fig3d))
message("Summary written to: ", file.path(results_dir, "E2_induced_transcription_SE_breaks_summary.csv"))
message("Figure written to: ", figures_dir)
