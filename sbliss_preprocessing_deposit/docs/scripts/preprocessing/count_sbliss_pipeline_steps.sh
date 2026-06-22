#!/usr/bin/env bash

set -euo pipefail

mkdir -p logs
out="logs/sbliss_pipeline_counts.tsv"
echo -e "step\tfile\tcount_type\tcount" > "$out"

count_fastq_reads() {
  local file="$1"
  echo $(( $(wc -l < "$file") / 4 ))
}

if compgen -G "processed_files/*.fastq" > /dev/null; then
  for file in processed_files/*.fastq; do
    echo -e "umi_extracted\t$file\treads\t$(count_fastq_reads "$file")" >> "$out"
  done
fi

if compgen -G "trimmed_output/*.fq" > /dev/null; then
  for file in trimmed_output/*.fq; do
    echo -e "trimmed\t$file\treads\t$(count_fastq_reads "$file")" >> "$out"
  done
fi

if compgen -G "mapped_output/*.sam" > /dev/null; then
  for file in mapped_output/*.sam; do
    echo -e "mapped_sam\t$file\talignments\t$(grep -v '^@' "$file" | wc -l)" >> "$out"
  done
fi

if compgen -G "bam_output/*.bam" > /dev/null; then
  for file in bam_output/*.bam; do
    echo -e "filtered_bam\t$file\talignments\t$(samtools view "$file" | wc -l)" >> "$out"
  done
fi

if compgen -G "deduplicated_output/*_deduplicated.bam" > /dev/null; then
  for file in deduplicated_output/*_deduplicated.bam; do
    echo -e "deduplicated_bam\t$file\tunique_reads\t$(samtools view "$file" | wc -l)" >> "$out"
  done
fi

if compgen -G "deduplicated_output/modified_*.bed" > /dev/null; then
  for file in deduplicated_output/modified_*.bed; do
    echo -e "break_site_bed\t$file\tbreak_sites\t$(wc -l < "$file")" >> "$out"
  done
fi

if compgen -G "bigwig_output/*.bedGraph" > /dev/null; then
  for file in bigwig_output/*.bedGraph; do
    echo -e "bedgraph\t$file\trows\t$(wc -l < "$file")" >> "$out"
  done
fi

if compgen -G "bigwig_output/*.bw" > /dev/null; then
  for file in bigwig_output/*.bw; do
    echo -e "bigwig\t$file\tfiles\t1" >> "$out"
  done
fi

if compgen -G "bigwigAverageOverBed_results/*.tab" > /dev/null; then
  for file in bigwigAverageOverBed_results/*.tab; do
    echo -e "gene_body_quantification\t$file\trows\t$(wc -l < "$file")" >> "$out"
  done
fi

echo "Pipeline count summary written to $out"
