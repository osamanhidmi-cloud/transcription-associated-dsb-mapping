#!/usr/bin/env python3
"""
UMI-aware deduplication for mapped sBLISS BAM files.

Reads are considered duplicates if they map to the same chromosome, start
position, and strand and carry the same UMI. The UMI is expected to be the final
field in the read name after splitting on '+'.
"""

import argparse
from pathlib import Path
import pysam


def remove_umi_duplicates(input_bam_path: Path, output_bam_path: Path):
    seen_umis = {}

    with pysam.AlignmentFile(input_bam_path, "rb") as input_bam, pysam.AlignmentFile(
        output_bam_path, "wb", template=input_bam
    ) as output_bam:
        for read in input_bam:
            if read.is_unmapped:
                continue

            read_name = read.query_name
            if "+" not in read_name:
                continue

            umi = read_name.split("+")[-1]
            position = (read.reference_name, read.reference_start, read.is_reverse)

            if position not in seen_umis:
                seen_umis[position] = set()

            if umi in seen_umis[position]:
                continue

            seen_umis[position].add(umi)
            output_bam.write(read)


def main():
    parser = argparse.ArgumentParser(description="Remove UMI duplicates from BAM files.")
    parser.add_argument("--input-dir", default="bam_output", help="Input BAM directory")
    parser.add_argument("--output-dir", default="deduplicated_output", help="Output BAM directory")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for bam_file in sorted(input_dir.glob("*.bam")):
        output_path = output_dir / bam_file.name.replace(".bam", "_deduplicated.bam")
        print(f"Processing {bam_file} -> {output_path}")
        remove_umi_duplicates(bam_file, output_path)


if __name__ == "__main__":
    main()
