# sBLISS preprocessing workflow

This document describes the raw-to-processed preprocessing workflow used to convert sBLISS FASTQ files into deduplicated break-site files, bigWig tracks, 10 kb bin-level signal tables, and gene-level signal matrices.

Raw FASTQ files are not included in this repository. They are available through GEO. The main downstream analysis scripts in this repository start from processed gene-level tables derived from the workflow below.

## Overview

The workflow consists of:

1. UMI/internal-index extraction and demultiplexing
2. Read-count check after UMI extraction
3. FastQC quality control
4. Adapter and quality trimming with Trim Galore
5. Alignment to the reference genome with HISAT2
6. SAM QNAME restoration to retain UMI/internal-index information
7. Filtering for uniquely mapped reads and conversion to sorted BAM
8. UMI-aware deduplication
9. Conversion from deduplicated BAM to BED
10. Conversion from read intervals to 1-bp break-site BED files
11. Generation of bedGraph and bigWig tracks
12. Quantification over 10 kb genomic bins
13. Quantification over gene bodies using `bigWigAverageOverBed`

## Required tools

The pipeline uses:

- Python 3
- FastQC
- Trim Galore
- HISAT2
- SAMtools
- BEDTools
- UCSC utilities: `bedGraphToBigWig`, `bigWigAverageOverBed`

Python scripts are provided in `scripts/preprocessing/`.

## Directory structure

The workflow expects or creates these directories:

```text
processed_files/
quality_reports/
trimmed_output/
mapped_output/
mapped_output_with_full_headers/
bam_output/
deduplicated_output/
bigwig_output/
bin_counts/
bigwigAverageOverBed_results/
logs/
```

## Step 1: Extract UMIs and internal indexes

The first 8 bases of each read are treated as the UMI. The next 8 bases are treated as the internal index. Reads are assigned to samples by matching the internal index sequence, allowing up to one mismatch by default.

Example:

```bash
python3 scripts/preprocessing/extract_umi.py \
  --output-dir processed_files \
  --indices CATCACGC GTCGTCGC ACGACCGC TGATGCGC CATCAATC GTCGTATC \
  --max-mismatches 1 \
  RPI1_S18_L008_R1_001.fastq.gz \
  RPI1_S18_L008_R2_001.fastq.gz
```

Count reads after UMI extraction:

```bash
echo "[Step 1] Read counts after UMI extraction:"
for file in processed_files/*.fastq; do
  echo "$file: $(($(wc -l < "$file") / 4)) reads"
done
```

## Step 2: Quality control

```bash
mkdir -p quality_reports
fastqc -o quality_reports processed_files/*.fastq
```

## Step 3: Adapter and quality trimming

```bash
mkdir -p trimmed_output

for file in processed_files/*.fastq; do
  trim_galore -q 20 --length 20 --output_dir trimmed_output "$file"
done
```

Count reads after trimming:

```bash
echo "[Step 3] Read counts after trimming:"
for file in trimmed_output/*.fq; do
  echo "$file: $(($(wc -l < "$file") / 4)) reads"
done
```

## Step 4: Build HISAT2 genome index

This step is only required once per reference genome.

```bash
hisat2-build genome.fa hisat2_index
```

## Step 5: Align reads using HISAT2

```bash
mkdir -p mapped_output
conda activate hisat2_env

for file in trimmed_output/*.fq; do
  basename=$(basename "$file" .fq)

  hisat2 -x hisat2_index -U "$file" \
    -S mapped_output/"$basename".sam \
    --no-spliced-alignment \
    --no-unal \
    --new-summary \
    --summary-file mapped_output/"$basename"_summary.txt
done
```

Count aligned reads:

```bash
echo "[Step 5] Aligned read counts, non-header SAM lines:"
for file in mapped_output/*.sam; do
  echo "$file: $(grep -v '^@' "$file" | wc -l) alignments"
done
```

## Step 6: Restore full SAM QNAMEs

Some aligners shorten read names in SAM output. This step restores the full read names from the trimmed FASTQ files, preserving the internal index and UMI added during UMI extraction.

```bash
python3 scripts/preprocessing/add_header_back.py \
  --mapped-dir mapped_output \
  --fastq-dir trimmed_output \
  --output-dir mapped_output_with_full_headers
```

## Step 7: Filter uniquely mapped reads and create sorted BAM files

Reads are retained if they have MAPQ ≥ 10 and the unique-alignment tag `NH:i:1`.

```bash
mkdir -p bam_output

for sam_file in mapped_output_with_full_headers/*.sam; do
  basename=$(basename "$sam_file" .sam)

  samtools view -h "$sam_file" | \
    awk 'BEGIN {OFS="\t"} /^@/ || ($5 >= 10 && $0 ~ /NH:i:1/) {print}' | \
    samtools view -b - | \
    samtools sort -o bam_output/"$basename".bam -

  samtools index bam_output/"$basename".bam
done
```

Count reads after filtering:

```bash
echo "[Step 7] Read counts in BAMs, MAPQ >= 10 and NH:i:1:"
for file in bam_output/*.bam; do
  echo "$file: $(samtools view "$file" | wc -l) alignments"
done
```

## Step 8: UMI-aware deduplication

```bash
python3 scripts/preprocessing/deduplicate.py \
  --input-dir bam_output \
  --output-dir deduplicated_output
```

Count unique reads after deduplication:

```bash
echo "[Step 8] Read counts after deduplication:"
for file in deduplicated_output/*.bam; do
  echo "$file: $(samtools view "$file" | wc -l) unique reads"
done
```

## Step 9: Index deduplicated BAM files

```bash
for file in deduplicated_output/*.bam; do
  samtools index "$file"
done
```

## Step 10: Convert deduplicated BAM files to BED

```bash
for file in deduplicated_output/*.bam; do
  basename=$(basename "$file" .bam)
  bedtools bamtobed -i "$file" > deduplicated_output/"$basename".bed
done
```

## Step 11: Create 1-bp break-site BED files

```bash
chmod +x scripts/preprocessing/create_break_bed.sh
scripts/preprocessing/create_break_bed.sh
```

## Step 12: Generate bedGraph files

```bash
mkdir -p bigwig_output

for bed in deduplicated_output/modified_*.bed; do
  base=$(basename "$bed" .bed)
  sort -k1,1 -k2,2n "$bed" > deduplicated_output/sorted_"$base".bed

  bedtools genomecov \
    -i deduplicated_output/sorted_"$base".bed \
    -g hg38.chrom.sizes.txt \
    -bg > bigwig_output/"$base".bedGraph
done
```

If chromosome naming differs between BAM/BED files and the chromosome-size file, create chromosome sizes from a BAM header:

```bash
samtools view -H example.bam \
  | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2"\t"$3}' \
  > genome_from_bam.chrom.sizes.txt
```

Then rerun `bedtools genomecov` with `genome_from_bam.chrom.sizes.txt`.

## Step 13: Convert bedGraph files to bigWig

```bash
for bg in bigwig_output/*.bedGraph; do
  base=$(basename "$bg" .bedGraph)

  bedGraphToBigWig \
    "$bg" \
    genome_from_bam.chrom.sizes.txt \
    bigwig_output/"$base".bw
done
```

## Step 14: Quantify signal over 10 kb genomic bins

Create 10 kb bins:

```bash
awk 'BEGIN{bin=0}
{
  chrom = $1
  size  = $2
  for (start = 0; start < size; start += 10000) {
    end = start + 10000
    if (end > size) end = size
    bin++
    printf("%s\t%d\t%d\tbin_%06d\n", chrom, start, end, bin)
  }
}' genome_from_bam.chrom.sizes.txt > bins10k.bed
```

Run `bigWigAverageOverBed`:

```bash
mkdir -p bin_counts

for bw in bigwig_output/*.bw; do
  base=$(basename "$bw" .bw)
  bigWigAverageOverBed "$bw" bins10k.bed bin_counts/"$base".tab
done
```

## Step 15: Count files and reads across preprocessing steps

```bash
chmod +x scripts/preprocessing/count_sbliss_pipeline_steps.sh
scripts/preprocessing/count_sbliss_pipeline_steps.sh
```

## Step 16: Quantify signal over gene bodies

```bash
mkdir -p bigwigAverageOverBed_results

for bw in bigwig_output/*.bw; do
  sample=$(basename "$bw" .bw)
  echo "Processing $sample"

  bigWigAverageOverBed \
    "$bw" \
    genes.sorted.bed \
    bigwigAverageOverBed_results/"${sample}.tab"
done
```

Create a gene-level signal matrix using the `sum` column from `bigWigAverageOverBed` output:

```bash
first_file=$(ls bigwigAverageOverBed_results/*.tab | head -n 1)

cut -f1 "$first_file" > gene_names.tmp
header="gene"

for file in bigwigAverageOverBed_results/*.tab; do
  sample=$(basename "$file" .tab)
  header="${header}\t${sample}"
  cut -f4 "$file" > "${sample}.sum.tmp"
done

{
  echo -e "$header"
  paste gene_names.tmp *.sum.tmp
} > genes_bigwig_sum_signal_matrix.tsv

rm gene_names.tmp *.sum.tmp
```

The resulting `genes_bigwig_sum_signal_matrix.tsv` can be used as input for downstream gene-level normalization, density calculation, and transcription-associated break analyses.

## Notes

This workflow documents the raw-to-processed sBLISS preprocessing strategy. Because raw sequencing files, reference indexes, and large intermediate files are not included in this GitHub repository, the main R analysis scripts start from processed gene-level tables.
