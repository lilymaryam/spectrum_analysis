import bte
import sys
import copy
import csv 
import random
import argparse
import random
from multiprocessing import Process
from collections import defaultdict
from scipy.stats import chi2_contingency

# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Process a phylogenetic tree to find splits, compute spectra, and get representative tips.")
    parser.add_argument("--input_tree", type=str, default="public-latest.all.masked.pb.gz", help="Input tree file (protobuf format)")
    parser.add_argument("--output_spectrum", type=str, default="spectra_output.tsv", help="Output TSV file for spectra")
    parser.add_argument("--min_chi", type=float, default=500, help="Minimum Chi-square value to accept a split")
    parser.add_argument("--min_mutations", type=int, default=500, help="Minimum number of mutations required for a split")
    parser.add_argument("--ntips", type=int, default=5, help="Number of tips to retrieve for each split")
    parser.add_argument("--bootstrap_splits", type=int, default=0, help="Number of bootstrap replicates to attempt in defining splits")
    parser.add_argument("--bootstrap_spectra", type=int, default=0, help="Number of bootstrap replicates to attempt in defining spectra")
    parser.add_argument("--nthreads", type=int, default=1, help="Number of threads for concurrent bootstrapping")
    parser.add_argument("--max_branch_length", type=int, default=100000, help="Maximum branch length to include in spectrum calculations")
    parser.add_argument("--calculate_max_chi", type=bool, default=False, help="Option to calculate maximum chi-square value based on tree size")
    return parser.parse_args()

### mutation positions 
def get_positions( node ) :
    positions = set() 
    for mutation in node.mutations:
        positions.add( int(mutation[1:-1]) )
    for child in node.children :
        positions.update( get_positions( child ) ) 
    return positions 

### create bootstrap weights by alignment position
def create_bootstrap( positions, n_samples=None ) :
    bootstrap_weights = defaultdict(int)  # Initialize a defaultdict to store counts
    if n_samples is None:
        n_samples = len(positions)  # Default to the size of the original set
    positions_list = list(positions)  # Convert set to list for indexing
    for _ in range(n_samples):
        selected_position = random.choice(positions_list)
        bootstrap_weights[selected_position] += 1
    return dict(bootstrap_weights)

def compute_mutation_spectrum(node, stop_nodes, spectrum_dict, weights=None, max_branch_length=100000):
    # Check if the spectrum for this node has already been computed
    if node in spectrum_dict:
        return spectrum_dict[node]
    
    #what is local spectrum?
    
    local_spectrum = defaultdict(int)
    #iterate throungh nodes of children
    for child in node.children:
        #what is a stop node?
        if any(child.id == stop_node.id for stop_node in stop_nodes):
            continue
        #recursive call to compute child spectrum
        child_spectrum = compute_mutation_spectrum(child, stop_nodes, spectrum_dict, weights, max_branch_length)
        #local spectrum is calculated as spectra of subtrees where current node is root
        for mutation_type, count in child_spectrum.items():
            local_spectrum[mutation_type] += count

    #why do we have max branch length ?
    #after getting lcal spectrum of children, we add current node to local spectrum
    if len(node.mutations) <= max_branch_length:

        for mutation in node.mutations:
            char1 = mutation[0]
            char2 = mutation[-1]
            mutation_type = char1 + char2
            pos = int(mutation[1:-1])
            
            # Assign weight 1 if weights is None; otherwise, check if pos is in weights dict and it gets that weigth or 0 otherwise
            weight = weights.get(pos, 0) if weights else 1
            local_spectrum[mutation_type] += weight
    
    #populate spectrum dict
    spectrum_dict[node] = local_spectrum
    return local_spectrum

def compute_spectrum_difference(spectrum1, spectrum2):
    difference_spectrum = defaultdict(int)
    all_keys = set(spectrum1.keys()).union(set(spectrum2.keys()))
    for key in all_keys:
        difference_spectrum[key] = abs(spectrum1.get(key, 0) - spectrum2.get(key, 0))
    return difference_spectrum

def subtract_spectra(spectrum1, spectrum2):
    remainder_spectrum = copy.deepcopy(spectrum1)
    for key, value in spectrum2.items():
        if key in remainder_spectrum:
            remainder_spectrum[key] -= value
    return remainder_spectrum

def normalize_spectrum(d):
    total = sum(d.values())
    if total == 0:
        raise ValueError("Cannot normalize because the total sum of values is 0.")
    normalized_d = {key: val / total for key, val in d.items()}
    return normalized_d

def get_spectra(finalized_splits, max_branch_length, weights=None):
    final_spectra = {} 
    for split_root in finalized_splits:
        print(f"Computing spectrum for subtree beginning at {split_root.id}", file=sys.stderr)
        spectrum_dict = {}
        final_spectra[split_root] = compute_mutation_spectrum(split_root, finalized_splits, spectrum_dict, weights, max_branch_length)
    return final_spectra

def write_spectra_to_tsv(spectra_dict, filename, ntips):
    all_keys = set()
    for spectrum in spectra_dict.values():
        all_keys.update(spectrum.keys())
    sorted_keys = sorted(all_keys)
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file, delimiter='\t')
        header = ["Node_ID"] + ["Total_Mutations"] + ["Number_Tips"] + ["Mutations:Tips"]+ sorted_keys + ["Exemplar tips"]
        writer.writerow(header)
        for node, spectrum in spectra_dict.items():
            tips = get_tips( spectra_dict.keys(), node )
            if ntips > 0 :
                normalized_spectrum = normalize_spectrum(spectrum)
                row = [node.id] + [sum(spectrum.values())] + [len(tips)] + [float(sum(spectrum.values()))/float(len(tips))] + [normalized_spectrum.get(key, 0) for key in sorted_keys] + [write_tips( tips, ntips )] 
                writer.writerow(row)
            else:
                normalized_spectrum = normalize_spectrum(spectrum)
                row = [node.id] + [sum(spectrum.values())] + [len(tips)] + [float(sum(spectrum.values()))/float(len(tips))] + [normalized_spectrum.get(key, 0) for key in sorted_keys]
                writer.writerow(row)
    print(f"Spectra written to {filename}", file=sys.stderr)

def get_tips(splits, node):
    tips = []
    def traverse(current_node):
        if current_node.is_leaf():
            tips.append(current_node.id)
        else:
            for child in current_node.children:
                if any(child.id == split.id for split in splits):
                    continue
                traverse(child)
    traverse(node)
    return tips

def write_tips(tips, ntips):
    if len(tips) == 0:
        return
    elif len(tips) < ntips:
        return ','.join(tips)
    else:
        return ','.join(random.sample(tips, ntips))

def find_splits(node, min_chi, min_mutations, max_branch_length, weights=None ):
    accepted_splits = set({node})
    finalized_splits = set()
    while len(accepted_splits) > len(finalized_splits):
        print(f"Starting iteration with {len(accepted_splits)-1} accepted splits and {len(finalized_splits)} finalized splits", file=sys.stderr)
        new_split = set()
        #iterate through current accepted splits
        for splitRoot in accepted_splits:
            #skip if already finalized
            if any(splitRoot.id == stop_node.id for stop_node in finalized_splits):
                print(f"Finalized split skipped:  {splitRoot.id}", file=sys.stderr)
                continue
            print(f"Computing spectrum for subtree beginning at {splitRoot.id}", file=sys.stderr)
            #spectrum dict will hold spectra for all subtrees with root [node] in tree
            spectrum_dict = {}
            #get spectrum for tree descending from node 
            #will also populate spectrum dict for all nodes in tree
            split_root_spectrum = compute_mutation_spectrum(splitRoot, accepted_splits, spectrum_dict, weights, max_branch_length)
            print(f"Computing distances between splits in (sub)tree {splitRoot.id}", file=sys.stderr)
            max_chi = 0
            max_chi_node = None
            #spectrum dict should have spectra for all subtrees in tree
            #node is no longer root of tree, but of subtree
            #iterate through all nodes in tree
            #this looks for subsequent possible splits within subtree rooted at splitRoot
            for node in spectrum_dict:
                #subtree from node needs to have enough mutations to be considered
                if sum(spectrum_dict[node].values()) < min_mutations:
                    continue
                #i beleive that this is the tree outside of the subtree rooted at node, may need to confirm
                #tree above splitroot ?
                #i feel like i dont get this. i feel like this should be spectrum of tree excluding subtree rooted at node
                #it seems like split_root_spectrum is spectrum of subtree rooted at splitRoot
                #i guess this doesnt matter a ton rn but i dont think this makes sense
                #####HELP
                ### this is confusing on 
                #splitroot sp
                if sum(split_root_spectrum.values()) - sum(spectrum_dict[node].values()) < min_mutations:
                    continue
                
                #note that meaning of splitroot_spectrum_list is spectrum of subtree rooted at node
                #meaning of above_spectrum_list is spectrum of tree rooted at splitRoot excluding subtree rooted at node
                splitroot_spectrum_list = []
                above_spectrum_list = [] 
                for mutation in ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"] :
                    if mutation in spectrum_dict[node].keys():
                        # Append the count from the node's spectrum, why?
                        splitroot_spectrum_list.append(spectrum_dict[node][mutation])
                    else:
                        # Append 0 if the mutation is missing from the spectrum
                        splitroot_spectrum_list.append(0)
                    
                    if mutation in split_root_spectrum.keys():
                        # Calculate the difference, but handle cases where mutation is missing from node
                        above_spectrum_list.append(split_root_spectrum[mutation] - spectrum_dict[node].get(mutation, 0))
                    else:
                        # Append 0 if the mutation is missing from the split_root_spectrum
                        above_spectrum_list.append(0)

                chi, p, dof, expected = chi2_contingency([splitroot_spectrum_list,above_spectrum_list])

                '''
                node_spectrum_difference = compute_spectrum_difference(
                    normalize_spectrum(spectrum_dict[node]),
                    normalize_spectrum(subtract_spectra(split_root_spectrum, spectrum_dict[node]))
                )
                chi = sum(node_spectrum_difference.values()) * min(
                    sum(spectrum_dict[splitRoot].values()) - sum(spectrum_dict[node].values()),
                    sum(spectrum_dict[node].values())
                )
                '''
                #check if chi is the largest seen so far
                #will only consider the split with the largest chi
                if chi > max_chi:
                    max_chi = chi
                    max_chi_node = node
            #after iterating through all nodes in tree, check if max chi is above threshold
            #if it is, add to new splits
            if max_chi > min_chi:
                
                if max_chi_node and max_chi_node not in accepted_splits:
                    new_split.add(max_chi_node)
                    print(f"New split found at {max_chi_node.id} with x2 {max_chi}", file=sys.stderr)
            else:
                finalized_splits.add(splitRoot)
                print(f"Finalized subtree rooted at {splitRoot.id}", file=sys.stderr)
        accepted_splits = accepted_splits.union(new_split)
        print(f"End of iteration: {len(new_split)} new splits added, {len(accepted_splits) -1} total accepted splits", file=sys.stderr)
    return finalized_splits

def bootstrap_replicate ( tree, replicate, min_chi, min_mutations, ntips, max_branch_length ) :
    print(f"Begining bootstrap no: {replicate}", file=sys.stderr)
    positions = get_positions( tree.root )
    bootstrap_weights = create_bootstrap( positions )
    finalized_splits_bootstrap = find_splits(tree.root, min_chi, min_mutations, max_branch_length, bootstrap_weights)
    bootstrap_spectra = get_spectra(finalized_splits_bootstrap, max_branch_length, bootstrap_weights)
    bootstrap_output_file = f"bootstrap_{replicate}_splits_output.tsv"
    write_spectra_to_tsv(bootstrap_spectra, bootstrap_output_file, ntips)

# Define the run_bootstrap function using explicit process creation
def run_bootstrap(tree, nbootstraps, nthreads, min_chi, min_mutations, max_branch_length):
    processes = []
    # Create and start a process for each bootstrap replicate
    for replicate in range(1, nbootstraps + 1):
        p = Process(target=bootstrap_replicate, args=(tree, replicate, min_chi, min_mutations, 0, max_branch_length))
        processes.append(p)
        p.start()
        # If we have reached the maximum number of threads, wait for them to finish
        if len(processes) >= nthreads:
            for p in processes:
                p.join()
            processes = []
    # Wait for any remaining processes to finish
    for p in processes:
        p.join()

    print(f"Bootstrap completed with {nbootstraps} replicates using {nthreads} threads.")

def bootstrap_spectrum_replicate( tree, replicate, splits, max_branch_lengths):
    print(f"Begining bootstrap no: {replicate}", file=sys.stderr)
    positions = get_positions( tree.root )
    bootstrap_weights = create_bootstrap( positions )
    bootstrap_spectra = get_spectra( splits, max_branch_lengths, bootstrap_weights)
    bootstrap_output_file = f"bootstrap_{replicate}_spectra_output.tsv"
    write_spectra_to_tsv(bootstrap_spectra, bootstrap_output_file, 0)

### ok, botostrap by spectrum
def run_bootstrap_spectra( tree, nbootstraps, nthreads, splits, max_branch_lengths ) :
    processes = []
    # Create and start a process for each bootstrap replicate
    for replicate in range(1, nbootstraps + 1):
        p = Process(target=bootstrap_spectrum_replicate, args=(tree, replicate, splits, max_branch_lengths))
        processes.append(p)
        p.start()
        # If we have reached the maximum number of threads, wait for them to finish
        if len(processes) >= nthreads:
            for p in processes:
                p.join()
            processes = []
    # Wait for any remaining processes to finish
    for p in processes:
        p.join()
    print(f"Bootstrap spectrum completed with {nbootstraps} replicates using {nthreads} threads.")

def main():

    ### read args and tree
    args = parse_args()
    tree = bte.MATree(args.input_tree)

    ### go through and do the real run without weighting mutations 
    finalized_splits = find_splits(tree.root, args.min_chi, args.min_mutations, args.max_branch_length )
    spectra = get_spectra(finalized_splits, args.max_branch_length )
    write_spectra_to_tsv(spectra, args.output_spectrum, args.ntips)

    ### get bootstrap splits if requested
    if ( args.bootstrap_splits > 0 ) :
        print(f"Bootstrapping splits with {args.bootstrap_splits} replicates using {args.nthreads} threads.", file=sys.stderr)
        run_bootstrap( tree, args.bootstrap_splits, args.nthreads, args.min_chi, args.min_mutations, args.max_branch_length )

    ### bootstrap spectrum requested:
    if ( args.bootstrap_spectra > 0 ) :
        print(f"Bootstrapping spectra with {args.bootstrap_spectra} replicates using {args.nthreads} threads.", file=sys.stderr)
        run_bootstrap_spectra( tree, args.bootstrap_spectra, args.nthreads, finalized_splits, args.max_branch_length )

if __name__ == "__main__":
    main()
