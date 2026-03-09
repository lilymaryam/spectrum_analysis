import numpy as np
from scipy.stats import chi2
import argparse
import bte

def bte_calculate(tree_file):
    tree = bte.MATree(tree_file)
    print(tree)
    tree.breadth_first_expansion()
    print('hi')
    internal_nodes = []
    mutations = {}
    #internal_nodes = [node.id for node in tree.breadth_first_expansion() if not node.is_leaf()]
    for node in tree.breadth_first_expansion():
        for mut in node.mutations:
            pos = mut[1:-1]
            #print(mut, pos)
            if pos not in mutations:
                mutations[pos] = 1
            else:
                mutations[pos] += 1
        if not node.is_leaf():
            internal_nodes.append(node.id)
    #print(internal_nodes)
    
    #muts = tree.count_mutation_types()
    print('muts', mutations)
    return internal_nodes, mutations

def mutation_cutoff(mutations, min_count):
    # We will use the 25th percentile (P25) to filter out the bottom 25% of least-frequent mutations.
    #percentile_threshold = 25
    mutation_counts = list(mutations.values())

    # Total number of unique mutation sites (This is your S before filtering)
    #s_initial = len(mutation_counts)
    s = 0 
    for count in mutation_counts:
        if count >= min_count:
            s += 1
    #print(f"Total unique mutation sites (S): {s}")
    return s
            


    # Use numpy.percentile to find the value that separates the bottom 25%
    #min_count_threshold = np.percentile(mutation_counts,  percentile_threshold)

    # Count how many sites remain after applying the new threshold (This is your final S)
    #s_final = sum(count >= min_count_threshold for count in mutation_counts)

    # --- 4. Output Results ---
    #print(f"Total unique mutation sites (S initial): {s_initial}")
    #print(f"Chosen percentile for filtering: P{percentile_threshold}")
    #print(f"Calculated --min_count threshold: {min_count_threshold}")
    #print("-" * 40)
    #print(f"Value to use for --min_count: {int(min_count_threshold)}")
    #print(f"Final number of sites tested (S final): {s_final}")


def get_mask_params(S, N, alpha_overall=0.05, DOF=1):
    # --- Input Parameters (You must change these to your actual counts!) ---
    #S_sites = 1000  # Estimated number of sites being tested
    #N_nodes = 50000 # Estimated number of internal nodes in the tree
    #alpha_overall = 0.05 # Desired overall False Positive rate (5%)
    #DOF = 1 # Degrees of Freedom for a 2x2 table

    # --- Bonferroni Calculation ---
    # 1. Calculate the total number of tests (T)
    T_total_tests = S * N
    print(T_total_tests)

    # 2. Calculate the Bonferroni-corrected p-value (alpha_prime)
    alpha_prime = alpha_overall / T_total_tests

    # 3. Find the Chi-squared critical value (X)
    # We use the Percent Point Function (ppf), which takes the area *to the left* (1 - alpha_prime)
    critical_chi2 = chi2.ppf(1 - alpha_prime, DOF)

    print(f"Estimated Sites (S): {S}")
    print(f"Estimated Nodes (N): {N}")
    print(f"Total Tests (T = S * N): {T_total_tests:,}")
    print(f"Bonferroni Corrected p-value (alpha'): {alpha_prime}")
    print(f"Required Chi-squared Critical Value (--mask_chi): {critical_chi2:.3f}")

def main():
    parser = argparse.ArgumentParser(description="Calculate mask parameters for spectrumSplits.")
    parser.add_argument("--tree", type=str, required=True, help="Estimated number of internal nodes in the tree")
    parser.add_argument("--min_count", type=int, required=False, default=20, help="Minimum count threshold for mutations")
    #parser.add_argument("--sites", type=int, required=True, help="Estimated number of sites being tested")
    args = parser.parse_args()
    int_nodes, muts = bte_calculate(args.tree)
    N = len(int_nodes)
    S = mutation_cutoff(muts, args.min_count)
    get_mask_params(S, N)
    #print('muts bitch', muts)
    #print('int nodes', int_nodes)
    #print(int_nodes) 
    #get_mask_params()

if __name__ == "__main__":
    main()
