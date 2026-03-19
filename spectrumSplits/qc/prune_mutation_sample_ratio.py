import bte
import sys
import argparse
import numpy as np

# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Process a phylogenetic tree for changepoint detection in mutation/descendant ratios.")
    parser.add_argument("--input_tree", type=str, required=True, help="Input tree file (protobuf format)")
    parser.add_argument("--output_tree", type=str, default="pruned_tree.pb.gz", help="Output tree file (protobuf format)")
    #need to fix this
    parser.add_argument("--threshold", type=float, default=0, help="Mutation:Leaf Ratio to prune. Note: If no value is set, or if value is set to 0, script will calculate best threshold. Recommended if unsure what mutation ratios are outliers.")
    parser.add_argument("--prune_list", type=str, default="to_prune.txt", help="File to write list of nodes to prune")
    return parser.parse_args()

# Compute the mutation-to-descendant ratio for each node and store it
def compute_descendants_mutations_ratio(node, mutation_ratio):
    if not node.children:  # Base case: If node is a tip (no children)
        mutation_ratio[node.id] = len(node.mutations) / 1  # 1 because it's a tip itself
        return 1, len(node.mutations), [node.id]  # Return the tip count and mutations

    total_tips = 0
    total_mutations = len(node.mutations)
    descendant_tips = []

    for child in node.children:  # Recursive case
        child_tips, child_mutations, child_descendant_tips = compute_descendants_mutations_ratio(child, mutation_ratio)
        total_tips += child_tips
        total_mutations += child_mutations
        descendant_tips.extend(child_descendant_tips)

    ratio = total_mutations / total_tips if total_tips > 0 else float('inf')
    mutation_ratio[node.id] = ratio

    return total_tips, total_mutations, descendant_tips

# Traverse the tree and detect changepoints based on mutation/descendant ratio changes
def detect_changepoints(node, mutation_ratio, threshold, changepoints, to_prune):
    for child in node.children:
        if mutation_ratio[child.id] >= threshold:
            changepoints.append((node.id, child.id, mutation_ratio[child.id], child))
            to_prune.add(child)
            continue
        detect_changepoints(child, mutation_ratio, threshold, changepoints, to_prune)

# Determine the threshold based on the overall distribution of ratios
def compute_threshold(mutation_ratios, k=3):
    ratios = np.array(list(mutation_ratios.values()))
    
    # Use Median Absolute Deviation (MAD) to set a robust threshold
    #NOTE: this doesn't work. code is here for posterity but mean + 3*stddev works much better
    #median = np.median(ratios)
    #abs_deviation = np.abs(ratios - median)
    #mad = np.median(abs_deviation)  
    #threshold = median + (k * mad)
    #print(f"Median Ratio: {median:.6f}")
    #print(f"MAD:          {mad:.6f}")
    #print(f"Threshold (k={k}): {threshold:.6f}")
    
    #This system is working to calculate mutation thresholds for each MAT
    #The original script for SARS-CoV-2 used mean + 2*stddev, which is better for larger datasets
    #But prunes too much for smaller datasets, 3*stddev is better for viral_usher_trees
    print(f"Mean mutation:descendant ratio: {np.mean(ratios)}")
    print(f"Stddev mutation:descendant ratio: {np.std(ratios)}")
    #print(f"Setting threshold to mean + 2*stddev: {np.mean(ratios) + np.std(ratios) * 2}")
    print(f"Setting threshold to mean + 3*stddev: {np.mean(ratios) + np.std(ratios) * 3}")
    return np.mean(ratios) + np.std(ratios) * 3
    #return threshold

# Function to prune marked nodes from the tree
def prune_tree(tree, to_prune):
    for node in to_prune:
        tree.remove_node(node.id)

# Helper function to get all descendant tips of a node
def get_descendant_tips(node):
    if not node.children:  # Base case: If node is a tip
        return [node.id]

    tips = []
    for child in node.children:
        tips.extend(get_descendant_tips(child))
    return tips

#Original version of script prunes internal nodes from tree, however this version outputs a list of samples to prune, which can be used to prune the tree in a separate step.
#This version is compatible with analyze.smk which will use matUtils to prune the tree based on the list of samples to prune. 
#This is because pruning internal nodes can potentially drastically change the tree and may also cause errors in BTE 

def main():
    args = parse_args()
    tree = bte.MATree(args.input_tree)

    mutation_ratio = {}
    compute_descendants_mutations_ratio(tree.root, mutation_ratio)

    #If threshold is set to 0 (default), script will compute a threshold 
    # based on the distribution of mutation:descendant ratios across the tree.
    if args.threshold == 0:
        args.threshold = compute_threshold(mutation_ratio)

    changepoints = []
    to_prune = set()
    detect_changepoints(tree.root, mutation_ratio, args.threshold, changepoints, to_prune)

    #printed output will be sent to logs in analyze.smk
    print("#Parent\tchild\ttips\tmutations:tips")
    #down stream analysis will be conducted on the file created 
    with open(args.prune_list, "w") as f:
        for parent_id, child_id, ratio, child_node in changepoints:
            descendant_tips = get_descendant_tips(child_node)
            tips_str = ",".join(descendant_tips) if descendant_tips else ""
            print(f"{parent_id}\t{child_id}\t{tips_str}\t{ratio}")
            #child_id is the internal node which is suggested for pruning. Pruning internal nodes instead of leaves
            #may cause errors. Pruning samples with matUtils will avoid errors
            #f.write(f"{child_id}\n")
            for tip in descendant_tips:
                f.write(f"{tip}\t{ratio}\t{args.threshold}\n")
    
if __name__ == "__main__":
    main()
