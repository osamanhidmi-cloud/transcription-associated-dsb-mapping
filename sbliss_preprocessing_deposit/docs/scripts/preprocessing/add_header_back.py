#!/usr/bin/env python3
"""
Restore full SAM QNAMEs from trimmed FASTQ files.

This is useful when alignment output contains shortened read names but the full
FASTQ read names contain UMI/internal-index information required for downstream
UMI-aware deduplication.
"""

import argparse
from pathlib import Path


def load_fastq_headers(fq_path: Path):
    """Map core QNAME to full QNAME from a FASTQ file."""
    qname_map = {}
    with open(fq_path, encoding="utf-8", errors="replace") as f:
        while True:
            line1 = f.readline().strip()
            if not line1:
                break
            f.readline()
            f.readline()
            f.readline()
            if line1.startswith("@"):
                full_qname = line1[1:]
                core_qname = full_qname.split("+")[0].split()[0]
                qname_map[core_qname] = full_qname
    return qname_map


def patch_sam_qnames(sam_path: Path, qname_map, output_path: Path):
    """Replace QNAMEs in SAM with full QNAMEs using a FASTQ-derived map."""
    with open(sam_path, encoding="utf-8", errors="replace") as fin, open(
        output_path, "w", encoding="utf-8"
    ) as fout:
        for line in fin:
            if line.startswith("@"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            core_qname = fields[0]
            if core_qname in qname_map:
                fields[0] = qname_map[core_qname]
            fout.write("\t".join(fields) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Restore full SAM QNAMEs using trimmed FASTQ files.")
    parser.add_argument("--mapped-dir", default="mapped_output", help="Directory containing SAM files")
    parser.add_argument("--fastq-dir", default="trimmed_output", help="Directory containing trimmed FASTQ/FQ files")
    parser.add_argument("--output-dir", default="mapped_output_with_full_headers", help="Output directory")
    args = parser.parse_args()

    mapped_dir = Path(args.mapped_dir)
    fastq_dir = Path(args.fastq_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for sam_file in mapped_dir.glob("*_trimmed.sam"):
        base = sam_file.stem
        fq_file = fastq_dir / f"{base}.fq"
        if not fq_file.exists():
            fq_file = fastq_dir / f"{base}.fastq"
        if not fq_file.exists():
            print(f"[Warning] No FASTQ found for: {base}")
            continue

        print(f"Processing: {sam_file.name} with {fq_file.name}")
        qname_map = load_fastq_headers(fq_file)
        output_sam = output_dir / sam_file.name
        patch_sam_qnames(sam_file, qname_map, output_sam)

    print("All SAM QNAMEs restored where matching FASTQ entries were found.")


if __name__ == "__main__":
    main()
