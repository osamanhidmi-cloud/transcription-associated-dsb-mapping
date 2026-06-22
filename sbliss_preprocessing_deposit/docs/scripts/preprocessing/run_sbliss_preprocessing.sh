#!/usr/bin/env bash

set -euo pipefail

# This script is a template for the sBLISS preprocessing workflow.
# Edit FASTQ_FILES, INTERNAL_INDICES, HISAT2_INDEX, and chromosome-size/BED paths before running.

FASTQ_FILES=(
  RPI1_S18_L008_R1_001.fastq.gz
  RPI1_S18_L008_R2_001.fastq.gz
  RPI2_S19_L008_R1_001.fastq.gz
  RPI2_S19_L008_R2_001.fastq.gz
  RPI3_S20_L008_R1_001.fastq.gz
  RPI3_S20_L008_R2_001.fastq.gz
  RPI4_S21_L008_R1_001.fastq.gz
  RPI4_S21_L008_R2_001.fastq.gz
  RPI5_S125_L008_R1_001.fastq.gz
  RPI5_S125_L008_R2_001.fastq.gz
  RPI6_S126_L008_R1_001.fastq.gz
  RPI6_S126_L008_R2_001.fastq.gz
)

INTERNAL_INDICES=(CATCACGC GTCGTCGC ACGACCGC TGATGCGC CATCAATC GTCGTATC)
HISAT2_INDEX="hisat2_index"
CHROM_SIZES="genome_from_bam.chrom.sizes.txt"
GENE_BED="genes.sorted.bed"

python3 scripts/preprocessing/extract_umi.py \
  --output-dir processed_files \
  --indices "${INTERNAL_INDICES[@]}" \
  --max-mismatches 1 \
  "${FASTQ_FILES[@]}"

mkdir -p quality_reports
fastqc -o quality_reports processed_files/*.fastq

mkdir -p trimmed_output
for file in processed_files/*.fastq; do
  trim_galore -q 20 --length 20 --output_dir trimmed_output "$file"
done

mkdir -p mapped_output
for file in trimmed_output/*.fq; do
  basename=$(basename "$file" .fq)
  hisat2 -x "$HISAT2_INDEX" -U "$file" \
    -S mapped_output/"$basename".sam \
    --no-spliced-alignment \
    --no-unal \
    --new-summary \
    --summary-file mapped_output/"$basename"_summary.txt
done

python3 scripts/preprocessing/add_header_back.py

mkdir -p bam_output
for sam_file in mapped_output_with_full_headers/*.sam; do
  basename=$(basename "$sam_file" .sam)
  samtools view -h "$sam_file" | \
    awk 'BEGIN {OFS="\t"} /^@/ || ($5 >= 10 && $0 ~ /NH:i:1/) {print}' | \
    samtools view -b - | \
    samtools sort -o bam_output/"$basename".bam -
  samtools index bam_output/"$basename".bam
done

python3 scripts/preprocessing/deduplicate.py

for file in deduplicated_output/*.bam; do
  samtools index "$file"
done

for file in deduplicated_output/*.bam; do
  basename=$(basename "$file" .bam)
  bedtools bamtobed -i "$file" > deduplicated_output/"$basename".bed
done

scripts/preprocessing/create_break_bed.sh

mkdir -p bigwig_output
for bed in deduplicated_output/modified_*.bed; do
  base=$(basename "$bed" .bed)
  sort -k1,1 -k2,2n "$bed" > deduplicated_output/sorted_"$base".bed
  bedtools genomecov -i deduplicated_output/sorted_"$base".bed -g "$CHROM_SIZES" -bg > bigwig_output/"$base".bedGraph
  bedGraphToBigWig bigwig_output/"$base".bedGraph "$CHROM_SIZES" bigwig_output/"$base".bw
done

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
}' "$CHROM_SIZES" > bins10k.bed

mkdir -p bin_counts
for bw in bigwig_output/*.bw; do
  base=$(basename "$bw" .bw)
  bigWigAverageOverBed "$bw" bins10k.bed bin_counts/"$base".tab
done

mkdir -p bigwigAverageOverBed_results
for bw in bigwig_output/*.bw; do
  sample=$(basename "$bw" .bw)
  bigWigAverageOverBed "$bw" "$GENE_BED" bigwigAverageOverBed_results/"${sample}.tab"
done

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

scripts/preprocessing/count_sbliss_pipeline_steps.sh
