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
    #print('r',root)
    #print('l', leaf.id)
    distance = 0
    current_node = leaf
    while current_node.id != root:
        #print(current_node.id)
        parent = current_node.parent
        #print('parent', parent.id)
        distance += len(current_node.mutations)
        #print(distance)
        current_node = parent
    return distance


def get_distances(leaf_id, bte_subtree, distance_reference):
    #for each leaf that is not the target leaf, get distance to target leaf
    #store distances to LCA avoid redundant calculations
    distance = {}
    #store list of distances in subtree for calculating average
    distances = []
    count = 0
    for leaf in bte_subtree.get_leaves():
        #leaf.id is the id of the sample compared to our target
        #leaf_id is the id of the target sample
        if leaf.id != leaf_id:
            LCA = bte_subtree.LCA([leaf_id, leaf.id])
            if leaf_id not in distance:
                leaf_id_dist = get_dist_to_root(LCA, bte_subtree.get_node(leaf_id))
                distance[leaf_id] = {LCA: leaf_id_dist}
            elif LCA not in distance[leaf_id]:
                leaf_id_dist = get_dist_to_root(LCA, bte_subtree.get_node(leaf_id))
                distance[leaf_id][LCA] = leaf_id_dist
            #distance from LCA to leaf_id node has been calculated before
            else:
                #print("time saved yay!")
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
            
            #print(leaf_id_dist, leaf_dist)
            distances.append(leaf_id_dist + leaf_dist)
            
            #distance = bte_subtree.get_distance(leaf_id, leaf.id)
            #yield (leaf.id, distance)
            #count +=1
        #if count == 1:
        #    break #for now
    return distances

def find_closest_neighbors(tree, portion=0.1, to_prune="to_prune.tsv"):
    #load full tree
    mat_tree = bte.MATree(tree)
    leaves = mat_tree.get_leaves()
    tree_size = len(leaves)
    
    #as this goes on, store distances to reduce redundant calculations
    dist_ref = {}
    
    #make a rescues directory for weird failures
    cmd = ["mkdir", "-p", "rescues"]
    subprocess.run(cmd, check=True)
    
    #calculate portion of tree size to extract 
    #LATER: have a max size to avoid huge extractions?
    portion_of_tree = int(tree_size * portion)
    #for each leaf, store the average distance to its closest X neighbors (X = portion_of_tree)
    #after all leaves are processed, find outliers in average distances 
    avg_distances = {}
    print(f"Tree size: {tree_size} portion: {portion} portion_of_tree: {portion_of_tree}")
    #parallelize this later?
    
    #count=0
    for leaf in leaves:
        #count +=1
        print(f"Processing leaf {count}/{tree_size}: {leaf.id}")
        #tempfile for analysis of each leaf
        with tempfile.TemporaryDirectory(prefix=f"{leaf.id}_") as tmpdir:
            #need a sample file for matUtils extract
            sample = os.path.join(tmpdir, f"{leaf.id}.txt")
            with open(sample, 'w') as f:
                f.write(f"{leaf.id}\n")
            #destination for extracted subtree around leaf
            out_path = os.path.join(tmpdir, f"{leaf.id}.pb.gz")
            #get subtree around leaf
            #will find distances of each sample to leaf in this subtree
            #thread restricted, may be worth modifying if used heavily
            #LATER: for each arg in comand create a variable and build cmd from those
            cmd = [
                "matUtils", "extract",
                "-i", tree,
                "-s", sample,
                "-Y", str(portion_of_tree),
                "-d", tmpdir,
                "-T", "1",
                "-o", f"{leaf.id}.pb.gz"
            ]
            try: 
                subprocess.run(cmd, check=True)
            #to catch when matUtils fails to extract
            #this happens rarely but needs to be handled
            except subprocess.CalledProcessError as e:
                #print('RESCUE')
                sample = os.path.join("./rescues/", f"{leaf.id}.txt")
                with open(sample, 'w') as f:
                    f.write(f"{leaf.id}\n")
                #thread restricted, may be worth modifying if used heavily
                cmd = [
                "matUtils", "extract",
                "-i", tree,
                "-s", sample,
                "-Y", str(portion_of_tree),
                "-d", "./rescues/",
                "-T", "1",
                "-o", f"{leaf.id}.pb.gz"
                ]
                print(f"Error running matUtils extract for leaf {leaf.id}: {e}")
                outpath = os.path.join("./rescues/", f"{leaf.id}.pb.gz")
                
            #print('here')
            #print(out_path)
            #print("tmpdir contents:", os.listdir(tmpdir))
            #make BTE subtree from extracted pb.gz
            bte_subtree = bte.MATree(out_path)
            #print('here2')
            #go find avg distance of this leaf to 10% of data 
            distances = get_distances(leaf.id, bte_subtree, dist_ref)
            avg_distance = sum(distances) / len(distances)
            avg_distances[leaf.id] = avg_distance
            #print(distances)
            #print('end')
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
    #this is a value i might need to toggle
    cutoff = 3
    mask = z > cutoff
    high_outliers = {k: v for k, v, m in zip(keys, values, mask) if m}

    #print("High outliers (z > 2):")
    with open(to_prune, 'w') as f:
        f.write(f"strain\tavg_distance\n")
        for sample in high_outliers:
            f.write(f"{sample}\t{high_outliers[sample]}\n")
    #for k, v in high_outliers.items():
    #    print(k, v, sep='\t')
    return high_outliers

def main():
    args = parse_args()
    find_closest_neighbors(args.input_tree, portion=args.portion_of_neighbors, to_prune=args.to_prune)

if __name__ == "__main__":
    main()
