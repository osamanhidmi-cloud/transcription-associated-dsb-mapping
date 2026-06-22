# ============================================================
# 01_score_transcription_stress_genes.R
#
# Purpose:
# Rank genes by transcription-stress score using the scoring method
# described in the paper:
#   - empirical P value from rank
#   - inverse normal transformation
#   - unweighted Liptak Z-score combination
#   - min-max scaling
#   - power transformation
#
# Genes not present among the top 1000 genes in all four parameters
# are excluded.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

input_dir <- "data/example"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_top_genes <- 1000

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

clean_gene_symbol <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    na_if("") |>
    toupper()
}

# Convert signal density to empirical P value based on rank.
# Higher signal = lower empirical P value.
add_empirical_p_and_z <- function(df, signal_col, prefix) {
  signal_col <- rlang::ensym(signal_col)
  
  df |>
    arrange(desc(!!signal_col)) |>
    mutate(
      rank = row_number(),
      n_genes = n(),
      
      # Higher signal gets smaller empirical P value.
      empirical_p = rank / (n_genes + 1),
      
      # Avoid infinite z-scores.
      empirical_p = pmin(pmax(empirical_p, 1 / (n_genes + 1)), n_genes / (n_genes + 1)),
      
      # Z_k = Phi^-1(1 - p_k)
      z_score = qnorm(1 - empirical_p)
    ) |>
    rename(
      "{prefix}_rank" := rank,
      "{prefix}_empirical_p" := empirical_p,
      "{prefix}_z" := z_score
    ) |>
    select(-n_genes)
}

min_max_scale <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  if (x_max == x_min) {
    return(rep(0, length(x)))
  }
  
  (x - x_min) / (x_max - x_min)
}

# ------------------------------------------------------------
# Read input tables
# ------------------------------------------------------------

break_density <- read_csv(
  file.path(input_dir, "break_density_example.csv"),
  show_col_types = FALSE
)

drip_density <- read_csv(
  file.path(input_dir, "drip_density_example.csv"),
  show_col_types = FALSE
)

top1_density <- read_csv(
  file.path(input_dir, "top1_density_example.csv"),
  show_col_types = FALSE
)

top1cc_density <- read_csv(
  file.path(input_dir, "top1cc_density_example.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Standardize signal tables
# ------------------------------------------------------------

# DSB signal from sBLISS.
# SE_mean represents the control condition:
# scramble siRNA + empty vector.
dsb_tbl <- break_density |>
  transmute(
    gene = clean_gene_symbol(symbol),
    dsb_density = SE_mean
  ) |>
  filter(!is.na(gene), !is.na(dsb_density)) |>
  group_by(gene) |>
  summarise(
    dsb_density = max(dsb_density, na.rm = TRUE),
    .groups = "drop"
  )

# R-loop signal from DRIP-seq.
rloop_tbl <- drip_density |>
  transmute(
    gene = clean_gene_symbol(symb),
    rloop_density = dens
  ) |>
  filter(!is.na(gene), !is.na(rloop_density)) |>
  group_by(gene) |>
  summarise(
    rloop_density = max(rloop_density, na.rm = TRUE),
    .groups = "drop"
  )

# TOP1 occupancy.
top1_tbl <- top1_density |>
  transmute(
    gene = clean_gene_symbol(symb),
    top1_density = dens
  ) |>
  filter(!is.na(gene), !is.na(top1_density)) |>
  group_by(gene) |>
  summarise(
    top1_density = max(top1_density, na.rm = TRUE),
    .groups = "drop"
  )

# TOP1cc signal.
# SE is the control signal column in this table.
top1cc_tbl <- top1cc_density |>
  transmute(
    gene = clean_gene_symbol(gene),
    top1cc_density = SE
  ) |>
  filter(!is.na(gene), !is.na(top1cc_density)) |>
  group_by(gene) |>
  summarise(
    top1cc_density = max(top1cc_density, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Add empirical P values and z-scores for each dataset
# ------------------------------------------------------------

dsb_scored <- dsb_tbl |>
  add_empirical_p_and_z(dsb_density, "dsb")

rloop_scored <- rloop_tbl |>
  add_empirical_p_and_z(rloop_density, "rloop")

top1_scored <- top1_tbl |>
  add_empirical_p_and_z(top1_density, "top1")

top1cc_scored <- top1cc_tbl |>
  add_empirical_p_and_z(top1cc_density, "top1cc")

# ------------------------------------------------------------
# Keep only genes among the top 1000 in all four parameters
# ------------------------------------------------------------

candidate_genes <- dsb_scored |>
  filter(dsb_rank <= n_top_genes) |>
  select(gene) |>
  inner_join(
    rloop_scored |> filter(rloop_rank <= n_top_genes) |> select(gene),
    by = "gene"
  ) |>
  inner_join(
    top1_scored |> filter(top1_rank <= n_top_genes) |> select(gene),
    by = "gene"
  ) |>
  inner_join(
    top1cc_scored |> filter(top1cc_rank <= n_top_genes) |> select(gene),
    by = "gene"
  ) |>
  distinct()

# ------------------------------------------------------------
# Combine scores using unweighted Liptak Z method
# ------------------------------------------------------------

transcription_stress_scores <- candidate_genes |>
  left_join(dsb_scored, by = "gene") |>
  left_join(rloop_scored, by = "gene") |>
  left_join(top1_scored, by = "gene") |>
  left_join(top1cc_scored, by = "gene") |>
  mutate(
    TSS_raw = (dsb_z + rloop_z + top1_z + top1cc_z) / sqrt(4),
    TSS_scaled = min_max_scale(TSS_raw),
    TSS_final = TSS_scaled^4
  ) |>
  arrange(desc(TSS_final)) |>
  mutate(TSS_rank = row_number())

# ------------------------------------------------------------
# Export results
# ------------------------------------------------------------

write_csv(
  transcription_stress_scores,
  file.path(output_dir, "transcription_stress_gene_scores.csv")
)

write_csv(
  transcription_stress_scores |> slice_head(n = 100),
  file.path(output_dir, "top_100_transcription_stress_genes.csv")
)

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

message("Analysis complete.")
message("Genes retained after top-1000 intersection: ", nrow(transcription_stress_scores))
message("Output written to: ", file.path(output_dir, "transcription_stress_gene_scores.csv"))
