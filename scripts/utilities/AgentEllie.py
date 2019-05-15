import subprocess
import argparse
from Bio import SeqIO
import pandas as pd

def iqtree(matrix, tree, threads):
    cmd = f"iqtree -s {matrix} -te {tree} -wsr -nt {threads} -m LG+F+R"
    print(cmd)
    subprocess.run(cmd, shell=True)

def parse_rates(rate_file):
    positions = []
    for line in open(rate_file):
        if "#" not in line and "Site" not in line:
            site, rate, _, _ = line.split('\t')
            positions.append((int(site)-1, float(rate)))
    sorted_positions = sorted(positions, key=lambda x: x[1], reverse=True)
    result = [i[0] for i in sorted_positions]
    return result

def main():
    # iqtree(args.matrix, args.tree, args.threads)
    matrix_dict = {}
    for record in SeqIO.parse(args.matrix, 'phylip'):
        matrix_dict[record.name] = pd.Series(list(record.seq))

    sorted_rates = parse_rates(args.matrix + ".rate")
    iter = 0
    for chunk in range(args.chunk, len(sorted_rates), args.chunk):
        with open(f'chunk{iter}', 'w') as res:
            res.write(f' {len(matrix_dict)} {len(sorted_rates) - chunk}\n')
            for name, seq in matrix_dict.items():
                res.write(f'{name} {"".join(seq[chunk:].values)}\n')
        iter += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='some description', usage="blabla")
    parser.add_argument('-m', '--matrix')
    parser.add_argument('-tr', '--tree')
    parser.add_argument('-c', '--chunk', type=int, default=3000)
    args = parser.parse_args()
    main()








