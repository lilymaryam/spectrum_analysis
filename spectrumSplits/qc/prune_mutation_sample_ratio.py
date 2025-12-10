import bte
import sys
import argparse
import numpy as np

# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Process a phylogenetic tree for changepoint detection in mutation/descendant ratios.")
    parser.add_argument("--input_tree", type=str, required=True, help="Input tree file (protobuf format)")
    parser.add_argument("--output_tree", type=str, default="pruned_tree.pb.gz", help="Output tree file (protobuf format)")
    parser.add_argument("--threshold", type=float, default=0, help="Mutation:Leaf Ratio to prune")
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
def compute_threshold(mutation_ratios):
    ratios = np.array(list(mutation_ratios.values()))
    print(f"Mean mutation:descendant ratio: {np.mean(ratios)}")
    print(f"Stddev mutation:descendant ratio: {np.std(ratios)}")
    print(f"Setting threshold to mean + 2*stddev: {np.mean(ratios) + np.std(ratios) * 2}")
    return np.mean(ratios) + np.std(ratios) * 2

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

def main():
    args = parse_args()
    tree = bte.MATree(args.input_tree)

    mutation_ratio = {}
    compute_descendants_mutations_ratio(tree.root, mutation_ratio)

    if args.threshold == 0:
        args.threshold = compute_threshold(mutation_ratio)

    changepoints = []
    to_prune = set()
    detect_changepoints(tree.root, mutation_ratio, args.threshold, changepoints, to_prune)

    print("#Parent\tchild\ttips\tmutations:tips")
    with open(args.prune_list, "w") as f:
        for parent_id, child_id, ratio, child_node in changepoints:
            descendant_tips = get_descendant_tips(child_node)
            tips_str = ",".join(descendant_tips) if descendant_tips else ""
            print(f"{parent_id}\t{child_id}\t{tips_str}\t{ratio}")
            #f.write(f"{child_id}\n")
            for tip in descendant_tips:
                f.write(f"{tip}\n")
    #print('made it')
    '''
    try:
        print('here')
        prune_tree(tree, to_prune)
        print('there')
        #tree.save_pb(args.output_tree)

        print('made it 1')
        
    except Exception as e:
        print(f"Error pruning tree: {e}")
    print('here')
    #print('done saving')
    try:
        tree.save_pb(args.output_tree)
    except Exception as e:
        print(f"Error saving tree: {e}")
    print('made it 2')
    '''
if __name__ == "__main__":
    main()
