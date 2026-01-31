import bte
import argparse
import subprocess
import tempfile
import os
import numpy as np

def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x!r} not a floating point literal")

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x!r} not in range [0.0, 1.0]")
    return x

def parse_args():
    parser = argparse.ArgumentParser(description="Prune samples from a BTE file based on a provided list.")
    parser.add_argument("--input_tree", "-i", required=True, type=str, help="Path to the input MAT file.")
    parser.add_argument("--portion_of_neighbors", "-n", type=restricted_float, default=0.1, help="Percentage of closest neighbors to find (0.0 to 1.0).")
    parser.add_argument("--to_prune", default="to_prune.tsv", type=str, help="Path to a text file containing sample names to keep, one per line.")
    #parser.add_argument("output_bte_file", type=str, help="Path to the output pruned BTE file.")
    return parser.parse_args()


def get_dist_to_root(root, leaf):
    print('r',root)
    print('l', leaf.id)
    distance = 0
    current_node = leaf
    while current_node.id != root:
        print(current_node.id)
        parent = current_node.parent
        print('parent', parent.id)
        distance += len(current_node.mutations)
        print(distance)
        current_node = parent
    return distance


def get_distances(leaf_id, bte_subtree):
    #for each leaf that is not the target leaf, get distance to target leaf
    print('here3')
    distance = {}
    distances = []
    count = 0
    for leaf in bte_subtree.get_leaves():
        #leaf.id is the id of the sample compared to our target
        #leaf_id is the id of the target sample
        if leaf.id != leaf_id:
            print('LeAF', leaf.id, 'leafid', leaf_id)

            LCA = bte_subtree.LCA([leaf_id, leaf.id])
            print('LCA',LCA)
            if leaf_id not in distance:
                leaf_id_dist = get_dist_to_root(LCA, bte_subtree.get_node(leaf_id))
                distance[leaf_id] = {LCA: leaf_id_dist}
            elif LCA not in distance[leaf_id]:
                leaf_id_dist = get_dist_to_root(LCA, bte_subtree.get_node(leaf_id))
                distance[leaf_id][LCA] = leaf_id_dist
            #distance from LCA to leaf_id node has been calculated before
            else:
                #time saved yay!
                leaf_id_dist = distance[leaf_id][LCA]
            
            if leaf.id not in distance:
                leaf_dist = get_dist_to_root(LCA, leaf)
                distance[leaf.id] = {LCA: leaf_dist}
            elif LCA not in distance[leaf.id]:
                leaf_dist = get_dist_to_root(LCA, leaf)
                distance[leaf.id][LCA] = leaf_dist
            #distance from LCA to leaf node has been calculated before
            else:
                #time saved yay!
                leaf_dist = distance[leaf.id][LCA]
            
            print(leaf_id_dist, leaf_dist)
            distances.append(leaf_id_dist + leaf_dist)
            
            #distance = bte_subtree.get_distance(leaf_id, leaf.id)
            #yield (leaf.id, distance)
            #count +=1
        #if count == 1:
        #    break #for now
    return distances

def find_closest_neighbors(tree, portion=0.01, to_prune="to_prune.tsv"):
    mat_tree = bte.MATree(tree)
    leaves = mat_tree.get_leaves()
    tree_size = len(leaves)
    portion_of_tree = int(tree_size * portion)
    avg_distances = {}
    print(f"Tree size: {tree_size} portion: {portion} portion_of_tree: {portion_of_tree}")
    #parallelize this later?
    for leaf in leaves:
        
        #tempfile for analysis of each leaf
        with tempfile.TemporaryDirectory(prefix=f"{leaf.id}_") as tmpdir:
            #need a sample file for matUtils extract
            print(tmpdir)
            sample = os.path.join(tmpdir, f"{leaf.id}.txt")
            with open(sample, 'w') as f:
                f.write(f"{leaf.id}\n")
            #destination for extracted subtree around leaf
            out_path = os.path.join(tmpdir, f"{leaf.id}.pb.gz")
            #get subtree around leaf
            #will find distances of each sample to leaf in this subtree
            cmd = [
                "matUtils", "extract",
                "-i", tree,
                "-s", sample,
                "-Y", str(portion_of_tree),
                "-d", tmpdir,
                "-o", f"{leaf.id}.pb.gz"
            ]

            subprocess.run(cmd, check=True)
            print('here')
            print(out_path)
            print("tmpdir contents:", os.listdir(tmpdir))
            bte_subtree = bte.MATree(out_path)
            print('here2')
            #go find avg distance of this leaf to 10% of data 
            distances = get_distances(leaf.id, bte_subtree)
            avg_distance = sum(distances) / len(distances)
            avg_distances[leaf.id] = avg_distance
            print(distances)
            print('end')
            #break
    #for leaf_id in avg_distances:
    #    print(f"{leaf_id}\t{avg_distances[leaf_id]}")
    '''
    arr = np.array(list(avg_distances.values()))
    z = (arr - arr.mean()) / arr.std()

    high_outliers = arr[z > 2] 
    print("High outliers (z > 2):", high_outliers)
    '''
    keys = list(avg_distances.keys())
    values = np.array(list(avg_distances.values()))

    z = (values - values.mean()) / values.std()

    mask = z > 2

    high_outliers = {k: v for k, v, m in zip(keys, values, mask) if m}

    print("High outliers (z > 2):")
    with open(to_prune, 'w') as f:
        for sample in high_outliers:
            f.write(f"{sample}\t{high_outliers[sample]}\n")
    #for k, v in high_outliers.items():
    #    print(k, v, sep='\t')
    return high_outliers




def main():
    args = parse_args()

    # Load the BTE file
    #find_closest_neighbors(args.input_tree, portion=args.portion_of_neighbors)

    high_outliers = find_closest_neighbors(args.input_tree, portion=args.portion_of_neighbors, to_prune=args.to_prune)
    
    # Read the samples to keep from the provided text file
    #with open(args.samples_to_keep, 'r') as f:
    #    samples_to_keep = {line.strip() for line in f if line.strip()}

    # Prune samples from the BTE data
    #pruned_bte_data = bte_data.prune_samples(samples_to_keep)

    # Save the pruned BTE data to the output file
    #pruned_bte_data.save(args.output_bte_file)

if __name__ == "__main__":
    main()
