Transcription-associated DNA break mapping analysis

This repository contains representative analysis workflows from my PhD work on endogenous transcription-associated DNA double-strand breaks (DSBs) in cancer cells.

The analyses focus on integrating gene-level DSB maps with transcription-associated genomic features, including R-loops, TOP1 occupancy, TOP1 cleavage complexes, super-enhancer-regulated genes, and estrogen-induced transcriptional activation.

The repository is intended to provide clean, readable, and runnable examples of the analysis logic. It is not a full raw sequencing data repository. Raw sequencing files and large intermediate files are not included.

Overview

The main analyses in this repository address three related questions:

1. Which genes show high transcription-associated DNA break potential when DSBs, R-loops, TOP1, and TOP1cc are integrated?
2. Are high transcription-stress genes enriched for super-enhancer-regulated genes and oncogenes?
3. Do super-enhancer-regulated genes show increased DNA break density across cancer and non-cancer cell lines?
4. Does estrogen-induced transcription increase DNA break density preferentially at ERSE/super-enhancer-associated genes?

Repository structure

transcription-associated-dsb-mapping/
  README.md
  R/
    01_score_transcription_stress_genes.R
    02_gene_set_enrichment.R
    03_example_plots.R
    04_SE_genes_cancer_vs_normal.R
    05_E2_induced_transcription_SE_breaks.R
  data/
    example/
      break_density_gene_level.csv
      rloop_density_gene_level.csv
      top1_density_gene_level.csv
      top1cc_density_gene_level.csv
      sedb_mcf7_genes.csv
      oncogenes.csv
      top_1000_expressed.csv
      figure1d_se_gene_break_enrichment.csv
      figure3d_E2_GROseq_SE_breaks.csv
  docs/
    sBLISS_preprocessing_workflow.md
  scripts/
    preprocessing/
      extract_umi.py
      add_header_back.py
      deduplicate.py
      create_break_bed.sh
      count_sbliss_pipeline_steps.sh
      run_sbliss_preprocessing.sh
  results/
  figures/

Data

The data/example/ directory contains processed gene-level input tables used by the R scripts.

Raw sequencing files are not included because of file size and data-sharing constraints. The raw sBLISS sequencing data are available through GEO. This repository starts from processed gene-level signal tables derived from the raw-to-processed sBLISS preprocessing workflow.

The included processed tables contain gene-level values such as:

* DSB density from sBLISS
* R-loop density from DRIP-seq
* TOP1 occupancy
* TOP1cc signal
* super-enhancer-regulated gene annotations
* oncogene annotations
* E2-induced GRO-seq fold-change groups
* E2-associated DSB densities before and after treatment

Raw sBLISS preprocessing

A raw-to-processed sBLISS preprocessing workflow is documented in:

docs/sBLISS_preprocessing_workflow.md

The documented workflow includes:

* UMI and internal-index extraction
* read quality control using FastQC
* adapter and quality trimming using Trim Galore
* alignment using HISAT2
* filtering for uniquely mapped reads
* UMI-aware deduplication
* conversion to break-site BED files
* generation of bedGraph and bigWig files
* quantification over 10 kb genomic bins
* quantification over gene bodies using bigWigAverageOverBed

The preprocessing scripts are provided in:

scripts/preprocessing/

These scripts document the raw-data processing logic, but the main runnable R analyses begin from processed gene-level tables.

R analysis scripts

01_score_transcription_stress_genes.R

This script integrates four gene-level transcription-stress-associated parameters:

1. DSB density
2. R-loop density
3. TOP1 occupancy
4. TOP1cc signal

For each dataset, genes are ranked by signal intensity. Ranks are converted to empirical P values and then to Z-scores. The four Z-scores are combined using an unweighted Liptak/Stouffer approach:

TSS_raw = (Z_DSB + Z_R-loop + Z_TOP1 + Z_TOP1cc) / sqrt(4)

The combined score is then min-max scaled and power transformed:

TSS_final = TSS_scaled^4

The output is a ranked table of transcription-stress genes.

Outputs:

results/transcription_stress_gene_scores.csv
results/top_100_transcription_stress_genes.csv

02_gene_set_enrichment.R

This script tests whether transcription-stress genes are enriched for selected gene sets using a hypergeometric test.

Gene sets include:

* super-enhancer-regulated genes
* oncogenes
* super-enhancer-regulated oncogenes
* highly expressed genes
* highly expressed super-enhancer-regulated genes

The gene universe is defined from the processed gene-level input tables, rather than only from the final transcription-stress gene list.

Output:

results/transcription_stress_gene_set_enrichment.csv

03_example_plots.R

This script generates simple example visualizations from the transcription-stress scoring and enrichment outputs.

The plotting script is intended as a basic output visualization step. It is not intended to exactly reproduce the final publication figure layout or styling.

Example outputs:

figures/top_transcription_stress_genes_barplot.pdf
figures/top_transcription_stress_genes_heatmap.pdf
figures/transcription_stress_gene_set_enrichment.pdf

04_SE_genes_cancer_vs_normal.R

This script analyzes DNA break enrichment at super-enhancer-regulated genes across cancer and non-cancer cell lines.

For each cell line, DNA break density at SE-regulated genes is represented relative to the median break density of non-SE genes from the same cell line.

This script uses a base R boxplot to stay close to the original analysis style.

Input:

data/example/figure1d_se_gene_break_enrichment.csv

Outputs:

results/SE_genes_cancer_vs_normal_summary.csv
figures/SE_genes_cancer_vs_normal_boxplot.pdf
figures/SE_genes_cancer_vs_normal_boxplot.png

05_E2_induced_transcription_SE_breaks.R

This script provides a cleaned implementation of the analysis logic underlying the Figure 3D-style comparison of E2-induced transcription and DNA break density.

Genes are grouped by their GRO-seq fold induction after E2 treatment. DNA break density is compared before and after E2 treatment for:

* non-ERSE / non-super-enhancer-associated genes
* ERSE / super-enhancer-associated genes

The script uses a corrected and cleaned implementation of the grouping and matched-sampling logic. Minor differences from the final published panel may occur because this repository emphasizes a transparent and corrected implementation of the analysis workflow.

Input:

data/example/figure3d_E2_GROseq_SE_breaks.csv

Outputs:

results/E2_induced_transcription_SE_breaks_summary.csv
figures/E2_induced_transcription_SE_breaks_boxplot.pdf
figures/E2_induced_transcription_SE_breaks_boxplot.png

How to run

Clone or download this repository, then open R or RStudio from the repository root.

Install required R packages if needed:

install.packages(c(
  "readr",
  "dplyr",
  "stringr",
  "purrr",
  "ggplot2",
  "pheatmap",
  "tibble"
))

Run the main analysis scripts:

source("R/01_score_transcription_stress_genes.R")
source("R/02_gene_set_enrichment.R")
source("R/03_example_plots.R")
source("R/04_SE_genes_cancer_vs_normal.R")
source("R/05_E2_induced_transcription_SE_breaks.R")

The scripts will create output files in:

results/
figures/

Notes on reproducibility

This repository is designed to make the gene-level analysis logic transparent and runnable.

The full raw-data workflow depends on large files that are not stored here, including:

* FASTQ files
* BAM files
* bedGraph files
* bigWig files
* genome FASTA files
* HISAT2 genome indexes
* chromosome-size files

The raw-to-processed strategy is documented in docs/sBLISS_preprocessing_workflow.md, while the R scripts use processed gene-level tables as inputs.

Software and external tools

The R analysis scripts use:

* R
* readr
* dplyr
* stringr
* purrr
* ggplot2
* pheatmap
* tibble

The preprocessing workflow uses external command-line tools including:

* FastQC
* Trim Galore
* HISAT2
* samtools
* bedtools
* UCSC tools such as bedGraphToBigWig and bigWigAverageOverBed
* Python
* pysam
* python-Levenshtein

Author

Osama Hidmi
PhD student, cancer research and genome instability
The Lautenberg Center for Immunology and Cancer Research
Hebrew University of Jerusalem
