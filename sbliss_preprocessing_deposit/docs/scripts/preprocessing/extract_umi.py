#!/usr/bin/env python3
"""
Extract UMIs and internal indexes from sBLISS FASTQ files.

The first 8 bases are treated as the UMI and the next 8 bases as the internal
index. Reads are assigned to output FASTQ files by internal-index sequence,
allowing a user-defined number of mismatches.
"""

import argparse
import gzip
from pathlib import Path


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def open_fastq(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def process_fastq_with_mismatches(input_files, output_dir, indices, max_mismatches=1):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    indices = list(indices)
    sample_files = {}
    match_counts = {idx: {"exact": 0, "mismatch": 0} for idx in indices}

    for input_file in input_files:
        input_file = Path(input_file)
        base_name = input_file.name.replace(".gz", "")
        for idx in indices:
            key = f"{base_name}_{idx}"
            sample_files[key] = open(output_dir / f"{key}.fastq", "w")

    for input_file in input_files:
        input_file = Path(input_file)
        print(f"Processing file: {input_file}")
        total_reads = 0
        matched_reads = 0
        base_name = input_file.name.replace(".gz", "")

        with open_fastq(input_file) as f:
            while True:
                line1 = f.readline().strip()
                if not line1:
                    break
                line2 = f.readline().strip()
                line3 = f.readline().strip()
                line4 = f.readline().strip()

                total_reads += 1

                if len(line2) < 16 or len(line4) < 16:
                    continue

                umi = line2[:8]
                internal_index = line2[8:16]

                new_line1 = f"{line1}+{internal_index}+{umi}"
                new_line2 = line2[16:]
                new_line4 = line4[16:]

                matched_index = None
                for idx in indices:
                    dist = hamming_distance(internal_index, idx)
                    if dist <= max_mismatches:
                        matched_index = idx
                        if dist == 0:
                            match_counts[idx]["exact"] += 1
                        else:
                            match_counts[idx]["mismatch"] += 1
                        break

                if matched_index:
                    key = f"{base_name}_{matched_index}"
                    matched_reads += 1
                    sample_files[key].write(
                        f"{new_line1}\n{new_line2}\n{line3}\n{new_line4}\n"
                    )

        print(
            f"File {input_file}: total reads processed = {total_reads}, "
            f"matched reads = {matched_reads}"
        )

    print("Match counts per index:")
    for idx, counts in match_counts.items():
        print(
            f"Index {idx}: exact matches = {counts['exact']}, "
            f"mismatched matches = {counts['mismatch']}"
        )

    for file_handle in sample_files.values():
        file_handle.close()


def main():
    parser = argparse.ArgumentParser(description="Extract UMI/internal-index sequences from FASTQ files.")
    parser.add_argument("fastq", nargs="+", help="Input FASTQ or FASTQ.gz files")
    parser.add_argument("--output-dir", default="processed_files", help="Output directory")
    parser.add_argument("--indices", nargs="+", required=True, help="Expected internal index sequences")
    parser.add_argument("--max-mismatches", type=int, default=1, help="Allowed mismatches to internal index")
    args = parser.parse_args()

    process_fastq_with_mismatches(
        input_files=args.fastq,
        output_dir=args.output_dir,
        indices=args.indices,
        max_mismatches=args.max_mismatches,
    )


if __name__ == "__main__":
    main()
