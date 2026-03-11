import os

configfile: "config.yaml"

rule all:
    input:
        expand("pruned/{virus}_final_pruned.jsonl.gz", virus=config["viruses"])
        #"pruned/Influenza_A_virus_A_California_07_2009_H1N1_4_final_pruned.jsonl.gz"
        #expand("pruned/{virus}_final_pruned.jsonl.gz", virus=config["good_viruses"])
        #
        #expand("reports/{virus}_split_report.txt", virus=config["good_viruses"])
        #expand("reports/{virus}_split_report.txt", virus=config["just_measles"]), expand("pruned/{virus}_final_pruned.jsonl.gz", virus=config["just_measles"])
        #expand("reports/{virus}_split_report.txt", virus=config["good_viruses"])
        #expand("pruned/{virus}_final_pruned.jsonl.gz", virus=config["good_viruses"])

rule prune:
    input:
        tree=os.path.join(config["data_dir"], "{virus}/optimized.pb.gz")
        #i probably need to change this to be whatever is the top scoring reroot
    output:
        prune_list="pruned/{virus}_to_prune.txt",
        #pruned_tree="pruned/{virus}_pruned.pb.gz"
    log:
        "logs/{virus}_prune.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p pruned
        python3 spectrumSplits/qc/prune_mutation_sample_ratio.py --input_tree {input.tree} --output_tree test --prune_list {output.prune_list} > {log} 2>&1
        """

rule matUtils:
    input:
        tree=os.path.join(config["data_dir"], "{virus}/optimized.pb.gz"), prunelist="pruned/{virus}_to_prune.txt"
    output:
        pruned_tree="pruned/{virus}_pruned.pb.gz", muts="data/{virus}_muts.txt"
    log:
        "logs/{virus}_matutils.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p redflags
        mkdir -p data
        matUtils extract -i {input.tree} -u pruned/{wildcards.virus}.samples > {log} 2>&1
        total=$(wc -l < pruned/{wildcards.virus}.samples)
        pruned=$(wc -l < {input.prunelist})
        echo "Total samples: $total" >> {log}
        echo "Samples to prune: $pruned" >> {log}
        echo "Pruning percentage: $(echo "scale=2; ($pruned / $total) * 100" | bc)% " >> {log}
        ratio=$(echo "$pruned / $total" | bc -l)
        if (( $(echo "$ratio > 0.1" | bc -l) )); then
            echo "More than 10% samples to be pruned, exiting." > redflags/{wildcards.virus}.txt
        fi
        if [ ! -s {input.prunelist} ]; then
            echo "Prunelist is empty, copying input tree to output." >> {log}
            cp {input.tree} {output.pruned_tree}
        else
            matUtils extract -i {input.tree} -s {input.prunelist} -p -o {output.pruned_tree} -T 1 >> {log} 2>&1
        fi
        matUtils summary -i {output.pruned_tree} -m data/{wildcards.virus}_muts.txt -T 1
        """
        



'''
rule look_for_weirdmuts:
    input:
        muts="data/{virus}_muts.txt"
    output:
        weird_muts="pruned/{virus}_weird_muts.txt"
    log:
        "logs/{virus}_weird_muts.log"
    shell:
        """
        python3 scripts/findweirdmuts.py -d {input.muts} -o {output.weird_muts} > {output.weird_muts} 
        """
'''

rule convert:
    input:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    output:
        final_tree="pruned/{virus}_final_pruned.jsonl.gz"
    log:
        "logs/{virus}_convert.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        slurm_extra="--export=ALL",
    shell:
        """
        usher_to_taxonium -i {input.pruned_tree} -t {wildcards.virus} -o {output.final_tree} > {log} 2>&1
        """

rule mask_muts:
    input:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    output:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz"
    log:
        "logs/{virus}_mask.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        slurm_extra="--export=ALL",
    shell:
        """
        python3 spectrumSplits/qc/mask_site_splits.py --input_tree {input.pruned_tree} --output_tree {output.masked_tree} > {log} 2>&1
        """

rule check_splits:
    input:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz"
    output:
        split_report="reports/{virus}_split_report.txt"
    log:
        "logs/{virus}_check_splits.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p reports
        python3 spectrumSplits/spectrumSplits/spectrumSplits.py --input_tree {input.masked_tree} --output_spectrum {output.split_report} > {log} 2>&1
        """
        