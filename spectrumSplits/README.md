# spectrumSplits
## Repository Contents
Spectrum splits is an algorithm for identifying changes in a method for subdividing a phylogeny into non-overlapping subtrees with distinct mutational spectra. This repository contains:  
  1. [Code and usage instructions for the main script.](https://github.com/russcd/spectrumSplits/tree/main/spectrumSplits)
  2. [Ancillary methods and usage instructions for automated QC of mutation annotated trees.](https://github.com/russcd/spectrumSplits/tree/main/qc)
  3. [Miscellaneous scripts for postprocessing and downstream analysis of results from the primary algorithm.](https://github.com/russcd/spectrumSplits/tree/main/misc)

## Dependencies
SpectrumSplits requires several libraries, all of which are available via conda and/or pip. The main dependency is the [Big-Tree-Explorer library](https://github.com/jmcbroome/bte-binder) for interacting with mutation-annotated-tree produced by [UShER](https://github.com/yatisht/usher)
