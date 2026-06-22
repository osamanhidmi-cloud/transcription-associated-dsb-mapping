# ============================================================
# 02_gene_set_enrichment.R
#
# Purpose:
# Test enrichment of gene sets among transcription-stress genes
# using the hypergeometric test.
#
# Important:
# The gene universe must be the genome-wide/background gene set,
# not only the final transcription-stress genes.
#
# Input:
#   results/transcription_stress_gene_scores.csv
#
# Files expected in data/example/:
#   break_density_gene_level.csv
#   rloop_density_gene_level.csv
#   top1_density_gene_level.csv
#   top1cc_density_gene_level.csv
#   sedb_mcf7_genes.csv
#   oncogenes.csv
#   top_1000_expressed.csv
#
# Output:
#   results/transcription_stress_gene_set_enrichment.csv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

input_dir <- "data/example"
results_dir <- "results"

# Since script 01 retained 64 genes, this will use all 64 if n_top > 64.
n_top <- 100
n_bottom <- 100

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

hypergeom_gene_set_test <- function(query_genes,
                                    gene_set,
                                    gene_universe,
                                    gene_set_name,
                                    direction = c("enrichment", "depletion")) {
  
  direction <- match.arg(direction)
  
  gene_universe <- unique(clean_gene_symbol(gene_universe))
  query_genes <- intersect(unique(clean_gene_symbol(query_genes)), gene_universe)
  gene_set <- intersect(unique(clean_gene_symbol(gene_set)), gene_universe)
  
  M <- length(gene_universe)                    # total background genes
  N <- length(query_genes)                      # selected genes
  m <- length(gene_set)                         # genes in gene set
  n <- M - m                                    # genes not in gene set
  k <- length(intersect(query_genes, gene_set)) # observed overlap
  
  expected <- (m / M) * N
  enrichment_ratio <- ifelse(expected == 0, NA_real_, k / expected)
  
  p_value <- if (direction == "enrichment") {
    phyper(k - 1, m, n, N, lower.tail = FALSE)
  } else {
    phyper(k, m, n, N, lower.tail = TRUE)
  }
  
  tibble(
    gene_set = gene_set_name,
    direction = direction,
    observed_overlap = k,
    expected_overlap = expected,
    enrichment_ratio = enrichment_ratio,
    p_value = p_value,
    query_size = N,
    gene_set_size = m,
    universe_size = M
  )
}

# ------------------------------------------------------------
# Read transcription-stress scores
# ------------------------------------------------------------

tss_scores <- read_csv(
  file.path(results_dir, "transcription_stress_gene_scores.csv"),
  show_col_types = FALSE
)

if (nrow(tss_scores) == 0) {
  stop("transcription_stress_gene_scores.csv is empty. Run script 01 first and check the input files.")
}

# ------------------------------------------------------------
# Read genome-wide processed signal tables to define background
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

# Background gene universe:
# all unique genes represented in the processed gene-level datasets.
# This replaces the hard-coded 26915 used in the exploratory script.
gene_universe <- Reduce(
  union,
  list(
    clean_gene_symbol(break_density$gene),
    clean_gene_symbol(rloop_density$gene),
    clean_gene_symbol(top1_density$gene),
    clean_gene_symbol(top1cc_density$gene)
  )
) |>
  unique() |>
  na.omit() |>
  as.character()

# ------------------------------------------------------------
# Define top and bottom transcription-stress gene groups
# ------------------------------------------------------------

n_top_use <- min(n_top, nrow(tss_scores))
n_bottom_use <- min(n_bottom, nrow(tss_scores))

top_tss_genes <- tss_scores |>
  arrange(desc(TSS_final)) |>
  slice_head(n = n_top_use) |>
  pull(gene)

bottom_tss_genes <- tss_scores |>
  arrange(TSS_final) |>
  slice_head(n = n_bottom_use) |>
  pull(gene)

# ------------------------------------------------------------
# Read gene-set files
# ------------------------------------------------------------

sedb_mcf7_genes <- read_csv(
  file.path(input_dir, "sedb_mcf7_genes.csv"),
  show_col_types = FALSE
)

oncogenes <- read_csv(
  file.path(input_dir, "oncogenes.csv"),
  show_col_types = FALSE
)

top_expressed <- read_csv(
  file.path(input_dir, "top_1000_expressed.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Standardize available gene sets
# ------------------------------------------------------------

se_genes <- sedb_mcf7_genes |>
  pull(gene) |>
  clean_gene_symbol()

oncogene_genes <- oncogenes |>
  pull(gene) |>
  clean_gene_symbol()

highly_expressed_genes <- top_expressed |>
  pull(gene) |>
  clean_gene_symbol()

se_oncogenes <- intersect(se_genes, oncogene_genes)
highly_expressed_se_genes <- intersect(highly_expressed_genes, se_genes)

# ------------------------------------------------------------
# Define gene sets to test
# ------------------------------------------------------------

gene_sets <- list(
  "SE-regulated genes" = se_genes,
  "Oncogenes" = oncogene_genes,
  "SE-regulated oncogenes" = se_oncogenes,
  "Highly expressed genes" = highly_expressed_genes,
  "Highly expressed SE-regulated genes" = highly_expressed_se_genes
)

# ------------------------------------------------------------
# Run enrichment/depletion tests
# ------------------------------------------------------------

top_enrichment <- imap_dfr(
  gene_sets,
  ~ hypergeom_gene_set_test(
    query_genes = top_tss_genes,
    gene_set = .x,
    gene_universe = gene_universe,
    gene_set_name = .y,
    direction = "enrichment"
  )
)

bottom_depletion <- imap_dfr(
  gene_sets,
  ~ hypergeom_gene_set_test(
    query_genes = bottom_tss_genes,
    gene_set = .x,
    gene_universe = gene_universe,
    gene_set_name = .y,
    direction = "depletion"
  )
)

enrichment_results <- bind_rows(top_enrichment, bottom_depletion) |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    minus_log10_p = -log10(p_value)
  ) |>
  arrange(direction, desc(enrichment_ratio))

# ------------------------------------------------------------
# Export results
# ------------------------------------------------------------

write_csv(
  enrichment_results,
  file.path(results_dir, "transcription_stress_gene_set_enrichment.csv")
)

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

message("Enrichment analysis complete.")
message("Background universe size: ", length(gene_universe))
message("Top query size: ", length(top_tss_genes))
message("Bottom query size: ", length(bottom_tss_genes))
message("Output written to: ", file.path(results_dir, "transcription_stress_gene_set_enrichment.csv"))
