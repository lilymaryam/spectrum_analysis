#!/usr/bin/env python3
"""Opportunity-normalize mononucleotide mutation spectra by genomic base composition.

Each mutation type X>Y can only occur at sites that are currently base X, so the
raw proportion of X>Y mutations is confounded by how much X the genome contains.
Dividing each type's proportion by the genomic frequency of its *source* base
removes that opportunity, giving a per-available-site rate. Rescaling the 12
corrected values back to sum to 1 yields spectra that are comparable across
genomes with different base compositions.
"""
import argparse
import sys
import pandas as pd

BASES = "ACGT"


def find_mutation_columns(cols):
    """Return columns named as an ordered pair of distinct bases (e.g. AC, GT)."""
    muts = []
    for c in cols:
        s = str(c).strip()
        if len(s) == 2 and s[0] in BASES and s[1] in BASES and s[0] != s[1]:
            muts.append(c)
    return muts


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-s", "--spectra", required=True,
                    help="Mutation spectra TSV (one row per sample, proportion per type).")
    ap.add_argument("-n", "--nucleotides", required=True,
                    help="Nucleotide count TSV with A,C,G,T columns.")
    ap.add_argument("-o", "--output", default="-",
                    help="Output TSV path (default: stdout).")
    ap.add_argument("--spectra-key", default="Sample",
                    help="Join column in the spectra file (default: Sample).")
    ap.add_argument("--nuc-key", default="Subdir",
                    help="Join column in the nucleotide file (default: Subdir).")
    ap.add_argument("--no-renormalize", action="store_true",
                    help="Emit raw opportunity-corrected values (proportion / source-base "
                         "frequency) without rescaling each row to sum to 1.")
    ap.add_argument("--keep-cols", default="Total_Mutations,Number_Tips",
                    help="Comma-separated spectra columns to carry through "
                         "(default: Total_Mutations,Number_Tips).")
    ap.add_argument("--report-freqs", action="store_true",
                    help="Also emit the genomic A/C/G/T frequencies used.")
    args = ap.parse_args()

    spec = pd.read_csv(args.spectra, sep="\t")
    nuc = pd.read_csv(args.nucleotides, sep="\t")

    for key, df, name in [(args.spectra_key, spec, "spectra"),
                          (args.nuc_key, nuc, "nucleotide")]:
        if key not in df.columns:
            sys.exit(f"Error: join column '{key}' not found in {name} file. "
                     f"Columns: {list(df.columns)}")

    mut_cols = find_mutation_columns(spec.columns)
    if len(mut_cols) != 12:
        sys.stderr.write(f"Warning: expected 12 mutation columns, found "
                         f"{len(mut_cols)}: {mut_cols}\n")
    if not mut_cols:
        sys.exit("Error: no mutation-type columns (e.g. AC, GT) found in spectra file.")

    for b in BASES:
        if b not in nuc.columns:
            sys.exit(f"Error: base column '{b}' not found in nucleotide file.")

    # Genomic base frequencies; N and Other are excluded from the denominator.
    denom = nuc[list(BASES)].sum(axis=1)
    freq = pd.DataFrame({b: nuc[b] / denom.where(denom > 0) for b in BASES})
    freq[args.nuc_key] = nuc[args.nuc_key].values

    merged = spec.merge(freq, left_on=args.spectra_key, right_on=args.nuc_key,
                        how="inner", suffixes=("", "_nucfreq"))

    unmatched = sorted(set(spec[args.spectra_key]) - set(merged[args.spectra_key]))
    if unmatched:
        preview = ", ".join(map(str, unmatched[:10])) + (" ..." if len(unmatched) > 10 else "")
        sys.stderr.write(f"{len(unmatched)} spectra sample(s) had no nucleotide match "
                         f"and were dropped: {preview}\n")

    # Opportunity normalization: divide each type by the frequency of its source base.
    norm = pd.DataFrame(index=merged.index)
    for m in mut_cols:
        src_freq = merged[m[0]]
        norm[m] = merged[m] / src_freq.where(src_freq > 0)  # NaN if source base absent

    if not args.no_renormalize:
        rowsum = norm.sum(axis=1)
        norm = norm.div(rowsum.where(rowsum > 0), axis=0)

    keep = [c for c in args.keep_cols.split(",") if c and c in merged.columns]
    pieces = [merged[[args.spectra_key] + keep]]
    if args.report_freqs:
        pieces.append(merged[list(BASES)].rename(columns={b: f"freq_{b}" for b in BASES}))
    pieces.append(norm)
    out = pd.concat(pieces, axis=1)

    dest = sys.stdout if args.output == "-" else args.output
    out.to_csv(dest, sep="\t", index=False, float_format="%.10g")
    if args.output != "-":
        sys.stderr.write(f"Wrote {len(out)} rows to {args.output}\n")


if __name__ == "__main__":
    main()
