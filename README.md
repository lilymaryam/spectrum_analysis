# spectrumSplits for viral usher trees 
This is a pipeline repository to analyze rerooted virus trees for mutation spectrum splits 

## Setup 
all dependencies are listed in envs/splits.yaml
spectrumSplits repo is expected to be downloaded or soft linked into the top dir of this repo.
I'm 99% sure usher is totally separate conda local build (after building you dont need the env active)

apparently usher for conda is fixed i need to check that at some point

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
