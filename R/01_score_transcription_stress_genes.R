# ============================================================
# 01_score_transcription_stress_genes.R
#
# Purpose:
# Rank genes by transcription-stress score using gene-level
# signal densities from:
#   1. endogenous DSBs / sBLISS
#   2. R-loops / DRIP-seq
#   3. TOP1 occupancy
#   4. TOP1cc signal
#
# Scoring method:
#   - Rank genes within each dataset
#   - Convert ranks to empirical P values
#   - Convert empirical P values to Z-scores
#   - Combine Z-scores using unweighted Liptak/Stouffer method
#   - Min-max scale the combined score
#   - Apply power transformation: TSS_final = TSS_scaled^4
#
# Input files expected in data/example/:
#   break_density_gene_level.csv
#   rloop_density_gene_level.csv
#   top1_density_gene_level.csv
#   top1cc_density_gene_level.csv
#   sedb_mcf7_genes.csv
#
# Output:
#   results/transcription_stress_gene_scores.csv
#   results/top_100_transcription_stress_genes.csv
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

# In the publication analysis, genes not among the top 1000
# in one or more parameters were excluded.
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

add_empirical_p_and_z <- function(df, signal_col, prefix) {
  signal_col <- rlang::ensym(signal_col)
  
  df |>
    arrange(desc(!!signal_col)) |>
    mutate(
      "{prefix}_rank" := row_number(),
      n_genes = n(),
      
      # Higher signal values correspond to lower empirical P values.
      "{prefix}_empirical_p" := .data[[paste0(prefix, "_rank")]] / (n_genes + 1),
      
      # Z_k = Phi^-1(1 - p_k)
      "{prefix}_z" := qnorm(1 - .data[[paste0(prefix, "_empirical_p")]])
    ) |>
    select(-n_genes)
}

min_max_scale <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  if (!is.finite(x_min) || !is.finite(x_max)) {
    stop("Cannot min-max scale: input contains no finite values.")
  }
  
  if (x_max == x_min) {
    return(rep(0, length(x)))
  }
  
  (x - x_min) / (x_max - x_min)
}

# ------------------------------------------------------------
# Read processed gene-level input tables
# ------------------------------------------------------------

break_density <- read_csv(
  file.path(input_dir, "break_density_gene_level.csv"),
  show_col_types = FALSE
)

rloop_density <- read_csv(
  file.path(input_dir, "rloop_density_gene_level.csv"),
  show_col_types = FALSE
)

top1_density <- read_csv(
  file.path(input_dir, "top1_density_gene_level.csv"),
  show_col_types = FALSE
)

top1cc_density <- read_csv(
  file.path(input_dir, "top1cc_density_gene_level.csv"),
  show_col_types = FALSE
)

sedb_mcf7_genes <- read_csv(
  file.path(input_dir, "sedb_mcf7_genes.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Standardize signal tables
# ------------------------------------------------------------

dsb_tbl <- break_density |>
  transmute(
    gene = clean_gene_symbol(gene),
    dsb_density = as.numeric(dsb_density)
  ) |>
  filter(!is.na(gene), !is.na(dsb_density)) |>
  group_by(gene) |>
  summarise(
    dsb_density = max(dsb_density, na.rm = TRUE),
    .groups = "drop"
  )

rloop_tbl <- rloop_density |>
  transmute(
    gene = clean_gene_symbol(gene),
    rloop_density = as.numeric(rloop_density)
  ) |>
  filter(!is.na(gene), !is.na(rloop_density)) |>
  group_by(gene) |>
  summarise(
    rloop_density = max(rloop_density, na.rm = TRUE),
    .groups = "drop"
  )

top1_tbl <- top1_density |>
  transmute(
    gene = clean_gene_symbol(gene),
    top1_density = as.numeric(top1_density)
  ) |>
  filter(!is.na(gene), !is.na(top1_density)) |>
  group_by(gene) |>
  summarise(
    top1_density = max(top1_density, na.rm = TRUE),
    .groups = "drop"
  )

top1cc_tbl <- top1cc_density |>
  transmute(
    gene = clean_gene_symbol(gene),
    top1cc_density = as.numeric(top1cc_density)
  ) |>
  filter(!is.na(gene), !is.na(top1cc_density)) |>
  group_by(gene) |>
  summarise(
    top1cc_density = max(top1cc_density, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Standardize SE-regulated gene list
# ------------------------------------------------------------

se_genes <- sedb_mcf7_genes |>
  transmute(gene = clean_gene_symbol(gene)) |>
  filter(!is.na(gene)) |>
  distinct() |>
  pull(gene)

# ------------------------------------------------------------
# Add empirical P values and Z-scores
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
# Keep genes among the top 1000 in all four parameters
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

if (nrow(candidate_genes) == 0) {
  stop(
    "No genes are shared among the top ",
    n_top_genes,
    " of all four parameters. Check whether gene identifiers match across datasets."
  )
}

# ------------------------------------------------------------
# Combine scores using unweighted Liptak/Stouffer method
# ------------------------------------------------------------

transcription_stress_scores <- candidate_genes |>
  left_join(dsb_scored, by = "gene") |>
  left_join(rloop_scored, by = "gene") |>
  left_join(top1_scored, by = "gene") |>
  left_join(top1cc_scored, by = "gene") |>
  mutate(
    TSS_raw = (dsb_z + rloop_z + top1_z + top1cc_z) / sqrt(4),
    TSS_scaled = min_max_scale(TSS_raw),
    TSS_final = TSS_scaled^4,
    SE_regulated = gene %in% se_genes
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
message("Genes retained after top-", n_top_genes, " intersection: ", nrow(transcription_stress_scores))
message("SE-regulated genes retained: ", sum(transcription_stress_scores$SE_regulated))
message("Output written to: ", file.path(output_dir, "transcription_stress_gene_scores.csv"))
