# spectrumSplits for viral usher trees 
This is a pipeline repository to analyze rerooted virus trees for mutation spectrum splits 

## Setup 
all dependencies are listed in envs/splits.yaml
spectrumSplits repo is expected to be downloaded or soft linked into the top dir of this repo.
I'm 99% sure usher is totally separate conda local build (after building you dont need the env active)

apparently usher for conda is fixed i need to check that at some point

# In this repo i have downloaded spectrumSplits (version? i dont think russ has updated it) and modified it for use with viral_usher_trees and my pipeline

## update 12/17/25
this is a fairly complex project consisting of 432 datasets and 2 (currently) pipelines. the first pipeline is outgrouping which is approximately 80% done with a few ongoing elements including 
- [ ] figure out if i can get more blast results for the 70 trees that dont have outgroups
- [ ] figure out how to classify good outgroups vs bad ones and 
- [ ] maybe retry outgrouping currently bad data
- [ ] figure out if the default ancestral node needs to be moved, when the outgroup is correct
- [ ] identify a good dataset to move forward with on spectrum splits (this will probably change multiple times but just to get a full pipeline)
- [ ] figure out how to stop qc prune mutation splits from overdoing it 
- [ ] figure out if mask sites is reasonable
- [ ] identify full tree spectra with spectrumsplits parameters that prevent iteration 
- [ ] run pca on full spectra


## Adapting spectrumSplits to viral_usher_trees
The original spectrumSplits algorithm relies on a chi2 cutoff of 500. This cutoff will not be appropriate for trees of smaller sizes (as all 432 viruses are). To adjust, 


## After running spectrumSplits
if using bootstraps, use process_bootstraps.py to analyze variation between (note i still dont know if bootstraps will be helpful )

update 1/26/26: 1000 bootstraps appear to support signal for measles. will check other trees as well.


## looking for more data.
human respirovirus 1 probably should be pruned before rerooting. will come back to this 

enterovirus d68 probably needs a more ancestral root but outgrouping looks alright. tree time for this runs into a similar issue 

## update 
I'm going to repeat spectra analysis on ../reroot_pipeline/viral_usher_trees/trees/Influenza_A_virus_A_California_07_2009_H1N1_4/timetree_rerooted.pb.gz. if it goes well i will try to generalize pipeline a bit 

## pipeline 
#prune weird samples based on errant mut rations (adjust mut ratios for different viruses?)

#general pipeline command 
`python3 spectrumSplits/qc/prune_mutation_sample_ratio.py --input_tree {input.tree} --output_tree test --prune_list {output.prune_list} > {log} 2>&1`
#command for flu
#note output tree doesnt currently exist. i disabled it in case it breaks.
#maybe reenable it i dont want to prune destructively
#im hoping these methods work well on good data, that is what im testing now
`python3 spectrumSplits/qc/prune_mutation_sample_ratio.py --input_tree ../reroot_pipeline/viral_usher_trees/trees/Influenza_A_virus_A_California_07_2009_H1N1_4/timetree_rerooted.pb.gz --output_tree Influenza_A_virus_A_California_07_2009_H1N1_4_pruned.pb.gz --prune_list Influenza_A_virus_A_California_07_2009_H1N1_4_prune_list`




