#!/usr/bin/env python3

import argparse
import glob
from collections import Counter
import os


def count_fasta(file_path):
    counts = Counter()

    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq = line.strip().upper()
            counts.update(seq)

    return counts


def main():
    parser = argparse.ArgumentParser(description="Count A/C/G/T per FASTA file")
    parser.add_argument(
        "--pattern",
        default="viral_usher_trees/**/*.fa",
        help="Glob pattern for FASTA files"
    )
    args = parser.parse_args()

    files = sorted(glob.glob(args.pattern, recursive=True))

    if not files:
        raise ValueError(f"No files found for pattern: {args.pattern}")

    # header
    print("Subdir\tSample\tA\tC\tG\tT\tN\tOther")

    for f in files:
        counts = count_fasta(f)

        A = counts.get("A", 0)
        C = counts.get("C", 0)
        G = counts.get("G", 0)
        T = counts.get("T", 0)
        N = counts.get("N", 0)
        other = sum(v for k, v in counts.items() if k not in "ACGTN")

        sample = os.path.basename(f)
        subdir = os.path.basename(os.path.dirname(f))  # <-- key addition

        print(f"{subdir}\t{sample}\t{A}\t{C}\t{G}\t{T}\t{N}\t{other}")


if __name__ == "__main__":
    main()
