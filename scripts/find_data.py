import sys 

file = sys.argv[1]

trees={}
with open(file, 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip('\n').split('\t')
        tree_name = parts[0]
        best_tree = parts[1]
        mid_r2 = parts[4]
        out_r2 = parts[10]
        treetime_r2 = parts[7]
        mid_slope = parts[5]
        out_slope = parts[11]
        treetime_slope = parts[8]
        data_available = [mid_r2, out_r2, treetime_r2, mid_slope, out_slope, treetime_slope]
        if best_tree == 'treetime':
            trees[tree_name] = {'method': 'treetime', 'r2': float(parts[7]), 'slope': float(parts[8])}
        elif best_tree == 'midpoint':
            trees[tree_name] = {'method': 'midpoint', 'r2': float(parts[4]), 'slope': float(parts[5])}
        elif best_tree == 'outgroup':
            trees[tree_name] = {'method': 'outgroup', 'r2': float(parts[10]), 'slope': float(parts[11])}
        else:
            trees[tree_name] = {'method': 'none', 'r2': 'NA', 'slope': 'NA'}
        
        if '' in data_available:
            print(data_available)
            trees[tree_name]['data_complete'] = False
        else:
            trees[tree_name]['data_complete'] = True


for tree in trees:
    #print(trees[tree])
    if trees[tree]['method'] != 'none':
        if trees[tree]['r2'] > 0.5:
            print(f"{tree}\t{trees[tree]['method']}\t{trees[tree]['r2']}\t{trees[tree]['slope']}")
    

        
