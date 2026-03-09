import argparse
import os

def check_muts_file(muts_file):
    positions = {}
    total = 0
    with open(muts_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            line = line.strip().split('\t')
            mut = line[0]
            position = mut[1:-1]
            count = int(line[1])
            total += count
            #print(mut,position, str(count))
            if position not in positions:
                positions[position] = count
            else:
                positions[position] += count
    #sorted_positions = dict(sorted(positions.values()))
    #print(sorted_positions)
    for p in positions:
        print(p, positions[p]/total)
    #    print(p, sorted_positions)

def main():
    parser = argparse.ArgumentParser(description="Find weird mutations in pruned trees.")
    parser.add_argument("--data", "-d", required=True, help="Path to file containing virus data.")
    parser.add_argument("--output", "-o", required=True, help="Output path for results.")
    args = parser.parse_args()

    check_muts_file(args.data)

if __name__ == "__main__":
    main()




