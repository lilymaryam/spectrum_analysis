# spectrumSplits analysis for viral_usher_trees

## Analysis
### Preprocessing 
#### prune errant samples 
can use `spectrumSplits/qc/prune_mutation_sample_ratio.py` as is used in the original SARS-CoV-2 analysis however, that dataset was much larger. 

`scripts/prune_samples.py` focuses on identifying samples that are outliers in terms of their distance to the closest 10% of samples in the dataset. running this script can be time consuming, however there are parameters that can be adjusted to improve runtime (i.e. size of dataset to compare to, anything else?)

#### mask aberrant mutations 
check for mutation positions that occur frequently and non uniformly across the tree. This code was written for SARS-CoV-2 and the default parameters are not going to work on the smaller viral datasets. For this reason, `--calculate_min_mutations` `--calculate_max_chi` should be added to command for running  mutation masking `python3 spectrumSplits/qc/mask_site_splits.py --input_tree {path to tree that has had errant samples pruned (previous step)}`

#### call spectrumSplits 
After pruning and masking, spectrum splits may be called `python3 spectrumSplits/spectrumSplits.py --input_tree {path to processed tree} --output_spectrum {path to spectrum} --calculate_min_chi`. `--calculate_min_chi` must be called to overwrite default cutoff values intended for the 8M sample SARS-CoV-2 tree.

##### Bootstraps
ill get to this

#### Visualizing spectra 
##### PCA
PCA can be visualized with `python3 spectrumSplits/misc/PCA.py -i {path to spectra tsv} -o {path to png of PCA figure}`

##### 
In order to visualize the output of the mutation spectra 2 scripts must be run. They are currently not in a pipeline but they will be added to one shortly.

first `python3 spectrumSplits/misc/annotate_nodes.py --spectrum_file {path to spectra tsv} --input_tree {path to preprocessed tree} --annotate_nodes_output_file {does the same thing as metadata?} --metadata_output {path to output metadata file to be used for visualization}`

then `python3 scripts/create_json.py --metadata_file {path to the metadata file created from annotate_nodes.py} --output_file {path to output json file to be used for visualization}`

**Side note: i get this error when i run the create_json.py script so that may need updating at some point. 
`/private/groups/corbettlab/lily/spectrumSplits/scripts/create_json.py:7: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead df = pd.read_csv(metadata_file, delim_whitespace=True)`**