import os
import random
import string

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

out_dir = config['out_dir']
in_format = config['in_format']
matrix = config['matrix']
tree = config['tree']
iqtree = config['iqtree']
rootfile = config['rootfile']
basename = config['basename']
rate_classes = config['rate_classes']

print(rate_classes)

def unique_name(keys):
    id_ = ''.join(random.choice(string.ascii_uppercase) for _ in range(10))
    if id_ not in keys:
        return id_
    else:
        unique_name(keys)


rule rename_seqfile:
    input:
        matrix
    output:
        f'{out_dir}/{basename}.renamed.phy',
        f'{out_dir}/{basename}.key.tsv'
    log:
        f'{out_dir}/logs/rename_seqfile/{basename}.log',
    run:
        pseudonames = {}
        records = []
        for record in SeqIO.parse(input[0], in_format):
            uname = unique_name(pseudonames)
            pseudonames[record.name] = uname
            records.append(SeqRecord(record.seq,
                                    id=uname,
                                    name='',
                                    description=''))
        SeqIO.write(records, output[0], 'phylip')
        
        with open(output[1], 'w') as outfile:
            index = 0
            for key, value in pseudonames.items():
                outfile.write(f'{key}\t{value}\t{index}\n')
                index +=1


rule rename_treefile:
    input:
        tree,
        f'{out_dir}/{basename}.key.tsv',
    output:
        f'{out_dir}/{basename}.renamed.tre',   
    run:
        with open(input[0], 'r') as in_treefile, open(input[1], 'r') as key_file, open(output[0], 'w') as out_treefile:
            tree_data = open(input[0]).read()
            for line in key_file:
                original_name, unique_name, _ = line.strip().split('\t')
                tree_data = tree_data.replace(original_name, unique_name)
            out_treefile.write(tree_data)


rule rename_iqtreefile:
    input:
        iqtree,
        f'{out_dir}/{basename}.key.tsv'
    output:
        f'{out_dir}/{basename}.renamed.iqtree'
    run:
        with open(input[0], 'r') as in_iqtreefile, open(input[1], 'r') as key_file, open(output[0], 'w') as out_iqtreefile:
            tree_data = open(input[0]).read()
            for line in key_file:
                original_name, unique_name, _ = line.strip().split('\t')
                tree_data = tree_data.replace(original_name, unique_name)
            out_iqtreefile.write(tree_data)


rule convert_rootfile:
    input:
        rootfile,
        f'{out_dir}/{basename}.key.tsv',
    output:
        f'{out_dir}/{basename}.rootfile.converted.txt'
    run:
        with open(input[1], 'r') as key_file:
            index_dict = {}
            for line in key_file:
                original_name, _, index = line.strip().split('\t')
                index_dict[original_name] = index

        with open(input[0], 'r') as in_rootfile, open(output[0], 'w') as out_rootfile:
            out_rootfile_frags = []
            for line in in_rootfile:
                taxon = line.strip()
                out_rootfile_frags.append(index_dict[taxon])

            out_rootfile.write(' '.join(out_rootfile_frags))
            out_rootfile.write('\n')


rule mammal:
    input:
        f'{out_dir}/{basename}.renamed.phy',
        f'{out_dir}/{basename}.renamed.tre',
    output:
        f'{out_dir}/{basename}.estimated-frequencies',
        f'{out_dir}/{basename}.esmodel.nex',
    params:
        rc=rate_classes
    conda:
        '../envs/mammal.yaml'
    shell: 
        '''
        mammal -s {input[0]} -t {input[1]} -c {params.rc} -l
        mv estimated-frequencies {output[0]}
        mv esmodel.nex {output[1]}
        '''


rule gfmix:
    input:
        matrix=f'{out_dir}/{basename}.renamed.phy',
        treefile=f'{out_dir}/{basename}.renamed.tre',
        iqtreefile=f'{out_dir}/{basename}.renamed.iqtree', 
        rootfile=f'{out_dir}/{basename}.rootfile.converted.txt'
    output:
        f'{out_dir}/{basename}.loglikelihood'
    params:
        aafreq='~/.gfmix/C20.aafreq.dat'
    conda:
        '../envs/gfmix.yaml'
    shell:
        '''
        gfmix -s {input.matrix} -t {input.treefile} -i {input.iqtreefile} -f {params.aafreq} -r {input.rootfile} > {output}
        '''