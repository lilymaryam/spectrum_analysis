import bte
import sys
import argparse
from scipy.stats import chi2_contingency
from collections import defaultdict
from multiprocessing import Process, Manager
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Process a phylogenetic tree to find splits, compute spectra, and mask mutations above/below nodes.")
    parser.add_argument("--input_tree", type=str, default="public-latest.all.masked.pb.gz", help="Input tree file (protobuf format)")
    parser.add_argument("--output_tree", type=str, default="masked_sites.pb.gz", help="Output tree file (protobuf format)")
    parser.add_argument("--min_total", type=int, default=500, help="Minimum mutation count to accept a split")
    parser.add_argument("--min_count", type=int, default=50, help="Minimum mutation count to check a site")
    parser.add_argument("--nthreads", type=int, default=100, help="Number of concurrent threads for processing")
    parser.add_argument("--mask_chi", type=float, default=5000, help="Minimum chi2 value for masking mutation below/above a node")
    return parser.parse_args()

def get_position_from_mutation(mutation):
    match = re.search(r'\d+', mutation)
    if match:
        return int(match.group(0))
    return None

def get_mutation_counts(node):
    mutation_counts = defaultdict(int)
    def traverse_tree(node):
        for mutation in node.mutations:
            pos = get_position_from_mutation(mutation)
            if pos is not None:
                mutation_counts[pos] += 1
        for child in node.children:
            traverse_tree(child)
    traverse_tree(node)
    return dict(mutation_counts)

def process_mutation(tree, position, count, total_mutations, args, mask_below_dict, mask_above_dict, chi_list):
    find_site_splits(position, count, total_mutations, tree.root, args, mask_below_dict, mask_above_dict, chi_list)

def run_in_process(tree, position, count, total_mutations, args, mask_below_dict, mask_above_dict, chi_list):
    p = Process(target=process_mutation, args=(tree, position, count, total_mutations, args, mask_below_dict, mask_above_dict, chi_list))
    p.start()
    return p

def find_site_splits(position, mutation_count, total_mutations, root, args, mask_below_dict, mask_above_dict, chi_list):
    mutation_memo = {}
    total_memo = {}
    max_chi = 0
    max_node = root
    mask_direction = 'below'

    def traverse_and_count(node):
        nonlocal max_chi, max_node, mask_direction

        if node.id in mutation_memo:
            return mutation_memo[node.id], total_memo[node.id]

        mutation_occurrences = sum(1 for mut in node.mutations if get_position_from_mutation(mut) == position)
        total_descendant_mutations = len(node.mutations)

        for child in node.children:
            child_mut_count, child_total = traverse_and_count(child)
            mutation_occurrences += child_mut_count
            total_descendant_mutations += child_total

        mutation_memo[node.id] = mutation_occurrences
        total_memo[node.id] = total_descendant_mutations

        snps_above = mutation_count - mutation_occurrences
        total_above = total_mutations - total_descendant_mutations

        rate_below = mutation_occurrences / total_descendant_mutations if total_descendant_mutations > 0 else 0
        rate_above = snps_above / total_above if total_above > 0 else 0

        if total_above > args.min_total and total_descendant_mutations > args.min_total:
            observed = [
                [total_above - snps_above, snps_above],
                [total_descendant_mutations - mutation_occurrences, mutation_occurrences]
            ]
            chi2, p, dof, expected = chi2_contingency(observed)

            if chi2 > max_chi:
                max_chi = chi2
                max_node = node
                mask_direction = 'below' if rate_below >= rate_above else 'above'

        return mutation_occurrences, total_descendant_mutations

    traverse_and_count(root)

    if max_chi > args.mask_chi:
        if mask_direction == 'below':
            current = mask_below_dict.get(max_node.id, [])
            current.append(position)
            mask_below_dict[max_node.id] = current
        else:
            current = mask_above_dict.get(max_node.id, [])
            current.append(position)
            mask_above_dict[max_node.id] = current

    chi_list.append((position, max_chi, max_node.id, mask_direction))

# ------------------ Masking logic ------------------

def mask_mutations(root, mask_below_dict, mask_above_dict):
    """Mask both below-node and above-node mutations."""

    # Helper to mask all descendants
    def mask_descendants(node, positions_to_mask):
        positions_to_mask = set(positions_to_mask)
        remaining_mutations = [m for m in node.mutations if get_position_from_mutation(m) not in positions_to_mask]
        node.update_mutations(remaining_mutations, update_branch_length=True)
        for child in node.children:
            mask_descendants(child, positions_to_mask)

    # Helper to mask everything except subtree of a node
    # Helper to mask everything except subtree of a node
    def mask_above(node, target_nodes):
        for target_id, positions in target_nodes.items():
            positions_set = set(positions)
            if not is_descendant(node, target_id) and node.id != target_id:
                remaining_mutations = [m for m in node.mutations if get_position_from_mutation(m) not in positions_set]
                node.update_mutations(remaining_mutations, update_branch_length=True)
        for child in node.children:  # keep your current traversal style
            mask_above(child, target_nodes)


    def is_descendant(node, ancestor_id):
        current = node
        while current.parent:
            if current.parent.id == ancestor_id:
                return True
            current = current.parent
        return False

    # Apply below-node masking
    for node_id, positions in mask_below_dict.items():
        node = find_node(root, node_id)
        if node:
            mask_descendants(node, positions)

    # Apply above-node masking
    mask_above(root, mask_above_dict)

def find_node(node, target_id):
    if node.id == target_id:
        return node
    for child in node.children:
        res = find_node(child, target_id)
        if res:
            return res
    return None

# ------------------ Main ------------------

def main():
    args = parse_args()
    tree = bte.MATree(args.input_tree)

    mutation_counts = get_mutation_counts(tree.root)
    print("Counting mutations", file=sys.stderr)
    total_mutations = sum(mutation_counts.values())

    mutation_counts = {k: v for k, v in mutation_counts.items() if v >= args.min_count}

    iteration = 1
    while len(mutation_counts) > 0:
        print(f"Finding splits. Iteration {iteration}", file=sys.stderr)

        manager = Manager()
        mask_below_dict = manager.dict()
        mask_above_dict = manager.dict()
        chi_list = manager.list()

        processes = []
        for pos, count in mutation_counts.items():
            print(f"\tPosition: {pos}\tOccurrences: {count}", file=sys.stderr)
            p = run_in_process(tree, pos, count, total_mutations, args, mask_below_dict, mask_above_dict, chi_list)
            processes.append(p)

            if len(processes) >= args.nthreads:
                for proc in processes:
                    proc.join()
                processes = []

        for proc in processes:
            proc.join()

        print(f"Mutations checked:", file=sys.stderr)
        for pos, chi, node_id, direction in sorted(chi_list, key=lambda x: -x[1]):
            print(f"{pos}\t{chi}\t{node_id}\t{direction}", file=sys.stderr)

        if args.mask_chi > 0 and (len(mask_below_dict) > 0 or len(mask_above_dict) > 0):
            print(f"Masking mutations: below={len(mask_below_dict)}, above={len(mask_above_dict)}", file=sys.stderr)
            mask_mutations(tree.root, mask_below_dict, mask_above_dict)

            # Recount mutations after masking
            mutation_counts = get_mutation_counts(tree.root)

            # Keep only positions that were masked for rechecking
            masked_positions = set(pos for positions in list(mask_below_dict.values()) + list(mask_above_dict.values()) for pos in positions)
            mutation_counts = {k: v for k, v in mutation_counts.items() if k in masked_positions}

            print(f"Sites to recheck: {len(mutation_counts)}", file=sys.stderr)
        else:
            mutation_counts = {}

        iteration += 1

    print(f"Saving tree to: {args.output_tree}", file=sys.stderr)
    tree.save_pb(args.output_tree)

if __name__ == "__main__":
    main()