#!/usr/bin/env bash

set -euo pipefail

mkdir -p logs

summary="logs/final_break_sites.tsv"
echo -e "sample\tinput_bed_rows\tbreak_bed_rows\ttotal_interval_length\tplus_breaks\tminus_breaks\tskipped_rows" > "$summary"

rm -f deduplicated_output/modified_*.bed
rm -f deduplicated_output/sorted_modified_*.bed

for bed in deduplicated_output/*_deduplicated.bed; do
  [ -e "$bed" ] || continue

  base=$(basename "$bed")

  case "$base" in
    modified_*|sorted_*) continue ;;
  esac

  out="deduplicated_output/modified_${base}"

  echo "Processing $bed"

  awk -v sample="${base%.bed}" -v summary="$summary" '
  BEGIN {
    OFS="\t"
    input_rows=0
    output_rows=0
    total_len=0
    plus=0
    minus=0
    skipped=0
  }

  {
    input_rows++
    chrom=$1
    start=$2
    end=$3

    # BED output may contain a read name with one space. In that case, the fields are:
    # $1 chrom, $2 start, $3 end, $4 read-name part 1, $5 read-name part 2,
    # $6 score, $7 strand.
    id=$4"_"$5
    score=1
    strand=$7

    if (strand == "+") {
      bstart=start
      bend=start+1
      plus++
    } else if (strand == "-") {
      bstart=end-1
      bend=end
      minus++
    } else {
      skipped++
      next
    }

    print chrom, bstart, bend, id, score, strand

    output_rows++
    total_len += bend - bstart
  }

  END {
    print sample, input_rows, output_rows, total_len, plus, minus, skipped >> summary

    if (output_rows == 0) {
      print "ERROR: no break sites were written for " sample > "/dev/stderr"
      exit 1
    }

    if (total_len != output_rows) {
      print "ERROR: intervals are not all 1 bp for " sample > "/dev/stderr"
      print "output_rows=" output_rows ", total_len=" total_len > "/dev/stderr"
      exit 1
    }

    if (skipped > 0) {
      print "WARNING: skipped " skipped " rows for " sample > "/dev/stderr"
    }
  }
  ' "$bed" > "$out"

  echo "Created $out"
done

echo "Summary written to $summary"
