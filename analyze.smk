import os

configfile: "config.yaml"

rule all:
    input:
        expand("pruned/{virus}_pruned.pb.gz", virus=config["viruses"][:5])

rule prune:
    input:
        tree=os.path.join(config["data_dir"], "{virus}/rerooted_outgroup_optimized.pb.gz")
    output:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    log:
        "logs/{virus}_prune.log"
    shell:
        """
        mkdir -p pruned
        python3 spectrumSplits/qc/prune_mutation_sample_ratio.py --input_tree {input.tree} --output_tree {output.pruned_tree} > {log} 2>&1
        """