# Automated qc tools
These scripts provide automated QC for reducing the impacts of potentially spurious mutations and samples. 

## prune_mutation_sample_ratio.py usage:
```usage: prune_mutation_sample_ratio.py [-h] [--input_tree INPUT_TREE] [--output_tree OUTPUT_TREE] [--threshold THRESHOLD]

Process a phylogenetic tree for changepoint detection in mutation/descendant ratios.

options:
  -h, --help            show this help message and exit
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --output_tree OUTPUT_TREE
                        Output tree file (protobuf format)
  --threshold THRESHOLD
                        Mutation:Leaf Ratio to prune
```

## mask_site_splits.py usage:
```usage: mask_site_splits.py [-h] [--input_tree INPUT_TREE] [--output_tree OUTPUT_TREE] [--min_total MIN_TOTAL] [--min_count MIN_COUNT] [--nthreads NTHREADS] [--mask_chi MASK_CHI]

Process a phylogenetic tree to mask highly imbalanced mutations.

options:
  -h, --help            show this help message and exit
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --output_tree OUTPUT_TREE
                        Output tree file (protobuf format)
  --min_total MIN_TOTAL
                        Minimum mutation count to accept a split
  --min_count MIN_COUNT
                        Minimum mutation count to accept a split
  --nthreads NTHREADS   Number of concurrent threads for processing
  --mask_chi MASK_CHI   Minimum chi2 value for masking mutation below a node (defaults to off)
```
