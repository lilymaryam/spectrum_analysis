# spectrumSplits analysis for viral_usher_trees
## Set Up 
This pipeline was run on a Linux HPC cluster. Snakemake and Slurm were used to distribute computing resources. 

usher was compiled with a conda local build v0.7.0. local build instructions explained here: https://usher-wiki.readthedocs.io/en/latest/Installation.html. 

conda environment was activated before running snakemake pipeline with slurm executor. `./envs/splits.yaml`

Need to download viral_usher_trees data. Git clone the github to a dir with sufficient space. Save the path to `viral_usher_trees/trees` in `data_dir` variable. Note that viral_usher_trees updates every month. Manually Git pull every month or automate data update. 

This codebase builds upon the code and analysis created by Russell Corbett-Detig in (spectrumSplits)[https://github.com/russcd/spectrumSplits]. I have copied relevant files from commit 74049f7. This is a clean copy. My spectrumSplits subdir contains scripts from this commit and modified for its own purposes. These files are no longer tracked to the original repo and are now a componenet of this repo and version controlled to this repository. 

## Analysis
These are the steps of the analysis explained. They are all connected through a snakemake pipeline `analyze.smk` which is used to run everything consecutively. 

### Preprocessing 
#### Prune errant samples 
Use `spectrumSplits/qc/prune_mutation_sample_ratio.py` as is used in the original SARS-CoV-2 analysis. This code is very similar to the original codebase with 2 main exceptions:

The first is that the threshold calculated for mutation/leaf ratio is 3 standard deviations above the mean instead of 2 (np.mean(ratios) + np.std(ratios) * 3) because these datasets are much smaller than the SARS-CoV-2 dataset.

The second is that the script no longer prunes the tree within the code using bte. This has been moved to a subsequent snakemake rule which uses matUtils extract for pruning to protect the pipeline from pruning edge cases in bte. 

#### Mask aberrant mutations 
After sample pruning, check for mutation positions that occur frequently and non uniformly across the tree. This code is adapted from `spectrumSplits/qc/mask_site_splits.py` in https://github.com/russcd/spectrumSplits (the original spectrumSplits database). It was written for SARS-CoV-2 and the default parameters are not going to work on the smaller viral datasets. 

For this reason, `--calculate_min_mutations` `--calculate_max_chi` are additional flags that have been added to the script and should be flagged in the command for running  mutation masking: `python3 spectrumSplits/qc/mask_site_splits.py --input_tree {path to tree that has had errant samples pruned (previous step)}` 

These 2 flags will change the minimum number of mutations needed to inspect a node (minimum will never be lower than 10) and chi threshold used to identify nodes where spectra changes. 

#### Call spectrumSplits 
After pruning and masking, spectrum splits may be called, however, the same issue persists where the original code makes assumptions about the dataset and hardcodes cutoffs. To adjust the algorithm to the datasets in viral_usher_trees  `--calculate_min_chi` must be called to overwrite default cutoff values intended for the 8M sample SARS-CoV-2 tree. `python3 spectrumSplits/spectrumSplits.py --input_tree {path to processed tree} --output_spectrum {path to spectrum} --calculate_min_chi`.

#### Identify single mutation spectra of tree from root node
For each viral tree, the mutation spectra of the root node (representing the tree as a whole) is calculating, restricting subsequent splits. The purpose of this is to get a sense of each dataset as a whole. These overall spectra are compared to each other in a PCA (`spectrumSplits/misc/PCA_single_split.py`). 

To do:
- [ ] add sars-c0v-2 pca to it. just take it from the old repo?
- [ ] filter bad datasets out 

##### Bootstraps
To examine the validity of the splits identified, bootstraps are highly recommended. I chose 1000. For some datasets, split identification will likely not be exact to the same node but similarity can be measured through jaccard similarity (spectra difference) and node distance. Bootstrap analysis is done with `/spectrumSplits/misc/process_bootstraps.py` which is highly similar to the original version found in https://github.com/russcd/spectrumSplits.


##### above this point is data generation. below this point are analyses which are still being workshopped. 


#### Visualizing spectra 
##### PCA
PCA can be visualized with `python3 spectrumSplits/misc/PCA.py -i {path to spectra tsv} -o {path to png of PCA figure}`

##### 
In order to visualize the output of the mutation spectra 2 scripts must be run. They are currently not in a pipeline but they will be added to one shortly.

first `python3 spectrumSplits/misc/annotate_nodes.py --spectrum_file {path to spectra tsv} --input_tree {path to preprocessed tree} --annotate_nodes_output_file {does the same thing as metadata?} --metadata_output {path to output metadata file to be used for visualization}`

then `python3 scripts/create_json.py --metadata_file {path to the metadata file created from annotate_nodes.py} --output_file {path to output json file to be used for visualization}`

**Side note: i get this error when i run the create_json.py script so that may need updating at some point. 
`/private/groups/corbettlab/lily/spectrumSplits/scripts/create_json.py:7: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead df = pd.read_csv(metadata_file, delim_whitespace=True)`**

after running these other commands, spectra changes on the tree can be visualized with `usher_to_taxonium --input {path to processed tree} --output {path to taxonium ready jsonl.gz} --metadata {path to metadata file to be used for visualization} -j {path to json config file from previous step} --columns SepctrumRoot,AC,AG,AT,CA,CG,CT,GA,GC,GT` the columnns are given to you exactly as they appear in the metadata . fix typo asap . 

## pipeline for analysis
To automate this process for viral_usher_trees analyze.smk runs the preprocessing, analysis, and visualization for all the tools. 

snakemake can be run on a server by delegating the number of cores with --cores {N cores}. For each rule, i have set resources that can be accessed by adding --slurm to the snakemake argument 
