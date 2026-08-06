# viralUSherSpectra

Tests whether a virus's mutation spectrum (rates of each of the 12 nucleotide
substitution types, e.g. A>C, C>T) is explained by its Baltimore classification.

Source trees: 
For this analysis, we used the March 2026 commit of `https://github.com/AngieHinrichs/viral_usher_trees` and the single split mutation spectra was calculated with `https://github.com/lilymaryam/spectrum_analysis` commit 1e90920. 

## Pipeline

Run in this order:

1. **`dedup.py`** — dedupes trees by sample/leaf overlap: loads tree files,
   ignores samples present in ~100% of the compared trees (reference/root
   artifacts), then keeps the larger tree whenever two trees' remaining
   samples overlap.
   ```
   dedup.py -i <tree_files...> > files2keep.txt
   ```
2. **`count_nts.py`** — counts A/C/G/T in each virus's rerooted reference fasta.
   → `nucleotide_counts.tsv`
3. *(external, not in this dir)* — raw per-type mutation tally from each
   deduped tree's branches (filtered to `files2keep.txt`) → raw 12-column
   spectrum per sample/segment. (THIS IS THE OUTPUT OF SPECTRUM_ANALYSIS)
4. **`normalize_spectra.py`** — opportunity-normalizes the raw spectrum:
   divides each mutation type's proportion by the genomic frequency of its
   source base (from `nucleotide_counts.tsv`), removing the base-composition
   confound.
   ```
   normalize_spectra.py -s <raw_spectrum.tsv> -n nucleotide_counts.tsv -o <normalized.tsv>
   ```
5. **`spectra_pca.py`** — collapses multi-segment viruses into one spectrum
   per virus, attaches Baltimore group (`baltimore_classification.tsv`), runs
   PCA, plots.
   ```
   spectra_pca.py -s <normalized.tsv> -b baltimore_classification.tsv \
       -o pca.png --html pca.html --dump-merged merged_dataset.tsv
   ```
6. **`permanova.py`** — distance-based PERMANOVA (default: Aitchison =
   CLR + Euclidean) testing how much spectrum variance Baltimore group
   explains, plus PERMDISP to flag dispersion confounds. Optionally add
   covariates from `virus_annotations.tsv`.
   ```
   permanova.py -i merged_dataset.tsv -a virus_annotations.tsv
   ```

## Files

| File | What |
|---|---|
| `dedup.py` | filters redundant/overlapping trees before analysis → `files2keep.txt` |
| `nucleotide_counts.tsv` | A/C/G/T counts per virus reference fasta |
| `baltimore_classification.tsv` | Baltimore group + genome type per virus (manually curated) |
| `virus_annotations.tsv` | host, enveloped y/n per virus |
| `merged_dataset.tsv` | final per-virus table: 12 spectrum cols + Baltimore group + metadata (output of `spectra_pca.py`) |
| `pca.png` / `pca.html` | PCA scatter, colored by Baltimore group |

## Caveats

- `merged_dataset.tsv` (162 viruses) was built from a **deduped tree
  snapshot** (273 trees, via `dedup.py`, from May). The live tree collection
  has since grown (446 trees) — rerun `dedup.py` against the current set
  before regenerating this analysis. Requires the `bte` and `rich` Python
  packages, and operates on tree protobufs
  (`viral_usher_trees/trees/*/optimized.pb.gz`), not the fasta files.
- `spectra_pca.py`'s PCA is **not** CLR-transformed (plain Euclidean on raw
  proportions); `permanova.py` **is** (Aitchison by default). The PCA plot is
  not drawn in the same geometry PERMANOVA tests in — don't read one as a
  visual proxy for the other.
- PERMANOVA treats each virus as independent. Closely related viruses within
  the same Baltimore group (e.g. Flaviviridae members all being group IV)
  aren't statistically independent — consider adding virus family as a term
  ahead of `Baltimore_group` in `--terms` before trusting a significant result.
