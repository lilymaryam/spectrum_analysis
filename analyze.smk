import os

configfile: "config.yaml"

rule all:
    input:
        expand("visuals/{virus}_split_visualization.jsonl.gz", virus=config["viruses"])
        #expand("reports/{virus}_bootstrap_report.tsv", virus=config["viruses"])
        #"reports/{virus}_bootstrap_report.txt
        #expand("reports/{virus}_multi_split_report.txt", virus=config["viruses"])
        #expand("reports/{virus}_single_split_report.txt", virus=config["viruses"])
        #expand("pruned/{virus}_final_pruned_masked.jsonl.gz", virus=config["viruses"])
        #"pruned/Measles_morbillivirus_final_pruned_masked.jsonl.gz"
        #expand("pruned/{virus}_pruned_masked.pb.gz", virus=config["viruses"])
        #do pruning only
        #expand("pruned/{virus}_final_pruned.jsonl.gz", virus=config["viruses"])

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
        #slurm_extra="--export=ALL",
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
    threads:
        config["threads"]
    log:
        "logs/{virus}_matutils.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
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
            matUtils extract -i {input.tree} -s {input.prunelist} -p -o {output.pruned_tree} -T {threads} >> {log} 2>&1
        fi
        matUtils summary -i {output.pruned_tree} -m data/{wildcards.virus}_muts.txt -T {threads}
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

#potentially can turn off this rule at some point 
rule convert_pruned:
    input:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    output:
        final_tree="pruned/{virus}_final_pruned.jsonl.gz"
    log:
        "logs/{virus}_convert_pruned.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        usher_to_taxonium -i {input.pruned_tree} -t {wildcards.virus} -o {output.final_tree} > {log} 2>&1
        """

rule mask_muts:
    input:
        pruned_tree="pruned/{virus}_pruned.pb.gz"
    output:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz"
    threads:
        config["threads"]
    log:
        "logs/{virus}_mask.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        python3 spectrumSplits/qc/mask_site_splits.py --input_tree {input.pruned_tree} --output_tree {output.masked_tree} --calculate_max_chi --calculate_min_mutations --nthreads {threads} > {log} 2>&1
        """

rule convert_masked:
    input:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz"
    output:
        final_tree="pruned/{virus}_final_pruned_masked.jsonl.gz"
    log:
        "logs/{virus}_convert_masked.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        usher_to_taxonium -i {input.masked_tree} -t {wildcards.virus} -o {output.final_tree} > {log} 2>&1
        """

rule check_single_split_spectra:
    input:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz"
    output:
        split_report="reports/{virus}_single_split_report.txt"
    threads:
        config["threads"]
    log:
        "logs/{virus}_check_single_split_spectra.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p bootstraps
        mkdir -p bootstraps/{wildcards.virus}
        mkdir -p reports
        python3 spectrumSplits/spectrumSplits.py --input_tree {input.masked_tree} --output_spectrum {output.split_report} --min_chi 100000 > {log} 2>&1
        """

rule single_split_pca:
    input:
        expand("reports/{virus}_single_split_report.txt", virus=config["viruses"]), metadata=os.path.join(config["data_dir"], "../tree_metadata.tsv")
    output:
        pca_plot="reports/single_split_pca.png",better_pca_plot="reports/single_split_pca_better.png"
    log:
        "logs/single_split_pca.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        python3 scripts/combine_single_splits.py --reports-dir ./reports --output ./reports/combined_single_splits.txt > logs/combine_single_splits.log 2>&1
        python3 spectrumSplits/misc/PCA_single_split.py --input reports/combined_single_splits.txt --output {output.better_pca_plot} --metadata {input.metadata} --results reports/single_split_pca_coordinates.txt > logs/single_split_pca_coordinates.log 2>&1
        """
        #python3 spectrumSplits/misc/PCA.py --input reports/combined_single_splits.txt --output {output.pca_plot} > {log} 2>&1

        
rule check_multi_split_spectra:
    input:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz"
    output:
        split_report="reports/{virus}_multi_split_report.txt"
    params:
        bootstrap_splits=config["bootstrap_replicates"]
    threads:
        config["threads"]
    log:
        "logs/{virus}_check_multi_split_spectra.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p bootstraps
        mkdir -p bootstraps/{wildcards.virus}
        mkdir -p reports
        python3 spectrumSplits/spectrumSplits.py --input_tree {input.masked_tree} --output_spectrum {output.split_report} --bootstrap_splits {params.bootstrap_splits} --bootstrap_dir bootstraps/{wildcards.virus} --calculate_min_chi > {log} 2>&1
        """

rule check_bootstraps:
    input:
        bootstrap_directory="bootstraps/{virus}", input_tree="pruned/{virus}_pruned_masked.pb.gz"
    output:
        bootstrap_report="reports/{virus}_bootstrap_report.tsv"
    params:
        bootstrap_splits=config["bootstrap_replicates"]
    threads:
        config["threads"]
    log:
        "logs/{virus}_check_bootstraps.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p reports
        python3 spectrumSplits/misc/process_bootstraps.py --bootstrap_directory {input.bootstrap_directory} --spectrum_file reports/{wildcards.virus}_multi_split_report.txt --output_file {output.bootstrap_report} --input_tree {input.input_tree} > {log} 2>&1
        """


rule visualize_splits:
    input:
        masked_tree="pruned/{virus}_pruned_masked.pb.gz", spectrum_file="reports/{virus}_multi_split_report.txt"
    output:
        metadata="visuals/{virus}_metadata.tsv", json="visuals/{virus}_metadata.json", split_viz="visuals/{virus}_split_visualization.jsonl.gz"
    log:
        "logs/{virus}_visualize_splits.log"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        #slurm_extra="--export=ALL",
    shell:
        """
        mkdir -p visuals
        python3 spectrumSplits/misc/annotate_nodes.py --spectrum_file {input.spectrum_file} --input_tree {input.masked_tree} --metadata_output {output.metadata}
        python3 scripts/create_json.py --metadata_file {output.metadata} --output_file {output.json}
        usher_to_taxonium --input {input.masked_tree} --output {output.split_viz}  --metadata {output.metadata}  -j {output.json} --columns SpectrumRoot,AC,AG,AT,CA,CG,CT,GA,GC,GT
        """
