import os

configfile: "config.yaml"

rule all:
    input:
        expand("reports/{virus}_split_report.txt", virus=config["good_viruses"])
        #expand("pruned/{virus}_final_pruned.jsonl.gz", virus=config["good_viruses"])

rule prune:
    input:
        tree=os.path.join(config["data_dir"], "{virus}/rerooted_no_outgroup_optimized.pb.gz")
    output:
        prune_list="pruned/{virus}_to_prune.txt",
        #pruned_tree="pruned/{virus}_pruned.pb.gz"
    log:
        "logs/{virus}_prune.log"
    shell:
        """
        mkdir -p pruned
        python3 spectrumSplits/qc/prune_mutation_sample_ratio.py --input_tree {input.tree} --output_tree test --prune_list {output.prune_list} > {log} 2>&1
        """

rule matUtils:
    input:
        tree=os.path.join(config["data_dir"], "{virus}/rerooted_no_outgroup_optimized.pb.gz"), prunelist="pruned/{virus}_to_prune.txt"
    output:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    log:
        "logs/{virus}_matutils.log"
    shell:
        """
        mkdir -p redflags
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
        echo "" >> {log}
        matUtils extract -i {input.tree} -s  {input.prunelist} -p -o {output.pruned_tree} >> {log} 2>&1
        """

rule convert:
    input:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    output:
        final_tree="pruned/{virus}_final_pruned.jsonl.gz"
    log:
        "logs/{virus}_convert.log"
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
    shell:
        """
        mkdir -p reports
        python3 spectrumSplits/spectrumSplits/spectrumSplits.py --input_tree {input.masked_tree} --output_spectrum {output.split_report} > {log} 2>&1
        """
        