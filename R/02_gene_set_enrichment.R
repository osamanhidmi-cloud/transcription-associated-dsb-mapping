# ============================================================
# 02_gene_set_enrichment.R
#
# Purpose:
# Test enrichment of gene sets among transcription-stress genes
# using the hypergeometric test.
#
# Important:
# The gene universe is defined from the full processed gene-level
# input tables, not from the final transcription-stress gene list.
#
# This script tests enrichment only. It does not perform depletion
# analysis, because the output of script 01 contains the retained
# high-scoring transcription-stress candidates, not a genome-wide
# ranked list with meaningful bottom-score genes.
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

hypergeom_enrichment_test <- function(query_genes,
                                      gene_set,
                                      gene_universe,
                                      gene_set_name) {
  
  gene_universe <- unique(clean_gene_symbol(gene_universe))
  query_genes <- intersect(unique(clean_gene_symbol(query_genes)), gene_universe)
  gene_set <- intersect(unique(clean_gene_symbol(gene_set)), gene_universe)
  
  M <- length(gene_universe)                       # background genes
  N <- length(query_genes)                         # selected TSS genes
  m <- length(gene_set)                            # genes in tested set
  n <- M - m                                       # genes not in tested set
  k <- length(intersect(query_genes, gene_set))    # observed overlap
  
  expected <- (m / M) * N
  enrichment_ratio <- ifelse(expected == 0, NA_real_, k / expected)
  
  # Upper-tail hypergeometric test:
  # probability of observing k or more overlaps by chance.
  p_value <- phyper(k - 1, m, n, N, lower.tail = FALSE)
  
  tibble(
    gene_set = gene_set_name,
    test = "enrichment_in_transcription_stress_genes",
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
  stop("transcription_stress_gene_scores.csv is empty. Run script 01 first.")
}

# Use all retained transcription-stress candidates from script 01.
tss_genes <- tss_scores |>
  arrange(desc(TSS_final)) |>
  pull(gene)

# ------------------------------------------------------------
# Read full processed signal tables to define gene universe
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

# Background universe:
# all genes represented in at least one processed gene-level table.
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
# Standardize gene sets
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

gene_sets <- list(
  "SE-regulated genes" = se_genes,
  "Oncogenes" = oncogene_genes,
  "SE-regulated oncogenes" = se_oncogenes,
  "Highly expressed genes" = highly_expressed_genes,
  "Highly expressed SE-regulated genes" = highly_expressed_se_genes
)

# ------------------------------------------------------------
# Run enrichment tests
# ------------------------------------------------------------

enrichment_results <- imap_dfr(
  gene_sets,
  ~ hypergeom_enrichment_test(
    query_genes = tss_genes,
    gene_set = .x,
    gene_universe = gene_universe,
    gene_set_name = .y
  )
) |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    minus_log10_p = -log10(p_value)
  ) |>
  arrange(desc(enrichment_ratio))

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
message("Transcription-stress query size: ", length(tss_genes))
message("Output written to: ", file.path(results_dir, "transcription_stress_gene_set_enrichment.csv"))
