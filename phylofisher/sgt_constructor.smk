# The main entry point of your workflow.
import os

import pandas as pd
from Bio import SeqIO

out_dir = config['out_dir']
in_dir = config['in_dir']
genes = config['genes'].split(',')
trees_only = config['trees_only']
no_trees = config['no_trees']

def get_output_files(wildcards):
    out_files=[]
    # I know the logic here is slightly confusing
    # It is this way because running alignment, other pre-tree processing, and trees are the default
    for gene in genes:
        if not no_trees:
            # If preprocessing and trees
            if not trees_only:
                with checkpoints.length_filtration.get(gene=gene).output[0].open() as infile:
                    line = infile.readline()
                    if line == '':
                        out_files.append(f'{out_dir}/length_filtration/bmge/{gene}.length_filtered')
                    else: 
                        out_files.append(f'{out_dir}/raxml/RAxML_bipartitions.{gene}.tre')
                        out_files.append(f'{out_dir}/raxml/RAxML_bootstrap.{gene}.tre')
                        out_files.append(f'{out_dir}/raxml/RAxML_bestTree.{gene}.tre')
                        out_files.append(f'{out_dir}/raxml/RAxML_info.{gene}.tre')
                        out_files.append(f'{out_dir}/raxml/RAxML_bipartitionsBranchLabels.{gene}.tre') 

            # If only trees
            else:
                out_files.append(f'{out_dir}/raxml/RAxML_bipartitions.{gene}.tre')
                out_files.append(f'{out_dir}/raxml/RAxML_bootstrap.{gene}.tre')
                out_files.append(f'{out_dir}/raxml/RAxML_bestTree.{gene}.tre')
                out_files.append(f'{out_dir}/raxml/RAxML_info.{gene}.tre')
                out_files.append(f'{out_dir}/raxml/RAxML_bipartitionsBranchLabels.{gene}.tre')              
        
        # If preprocessing but no trees
        else: 
            with checkpoints.length_filtration.get(gene=gene).output[0].open() as infile:
                    line = infile.readline()
                    if line == '':
                        out_files.append(f'{out_dir}/length_filtration/bmge/{gene}.length_filtered')
                    else: 
                        out_files.append(f'{out_dir}/trimal/{gene}.final')
            
    
    return out_files

rule all:
        input: get_output_files

# if not trees_only:
rule rm_star_gaps:
    input:
        f'{in_dir}/{{gene}}.fas'
    output:
        f'{out_dir}/prequal/{{gene}}.aa'
    run:
        with open(output[0], 'w') as res:
            for record in SeqIO.parse(input[0], 'fasta'):
                res.write(f'>{record.name}\n{str(record.seq).replace("-", "").replace("*", "")}\n')

rule prequal:
    input:
        f'{out_dir}/prequal/{{gene}}.aa'
    output:
        f'{out_dir}/prequal/{{gene}}.aa.filtered'
    log:
        f'{out_dir}/logs/prequal/{{gene}}.log'
    shell:
        'prequal {input} >{log} 2>{log}'

rule length_filter_mafft:
    input:
        f'{out_dir}/prequal/{{gene}}.aa.filtered'
    output:
        f'{out_dir}/length_filtration/mafft/{{gene}}.aln'
    log:
        f'{out_dir}/logs/length_filter_mafft/{{gene}}.log'
    shell:
        'mafft --thread 1 --globalpair --maxiterate 1000 --unalignlevel 0.6 {input} >{output} 2>{log}'

rule length_filter_divvier:
    input:
        f'{out_dir}/length_filtration/mafft/{{gene}}.aln'
    output:
        f'{out_dir}/length_filtration/divvier/{{gene}}.aln.partial.fas',
        f'{out_dir}/length_filtration/divvier/{{gene}}.aln.PP'
    log:
        f'{out_dir}/logs/length_filter_divvier/{{gene}}.log'
    shell:
        f'''
        divvier -mincol 4 -partial {{input}} >{{log}} 2>{{log}}

        mv {out_dir}/length_filtration/mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/length_filtration/divvier >{{log}} 2>{{log}}
        mv {out_dir}/length_filtration/mafft/{{wildcards.gene}}.aln.PP {out_dir}/length_filtration/divvier >{{log}} 2>{{log}}
        '''

rule x_to_dash:
    input:
        f'{out_dir}/length_filtration/divvier/{{gene}}.aln.partial.fas'
    output:
        f'{out_dir}/length_filtration/bmge/{{gene}}.pre_bmge'
    log: 
        f'{out_dir}/logs/x_to_dash/{{gene}}.log'
    run:
        with open(output[0], 'w') as res:
            for record in SeqIO.parse(input[0], 'fasta'):
                res.write(f'>{record.name}\n{str(record.seq).replace("X", "-")}\n')

rule length_filter_bmge:
    input:
        f'{out_dir}/length_filtration/bmge/{{gene}}.pre_bmge'
    output:
        f'{out_dir}/length_filtration/bmge/{{gene}}.bmge'
    log:
        f'{out_dir}/logs/length_filter_bmge/{{gene}}.log'
    shell:
        'bmge -t AA -g 0.3 -i {input} -of {output} >{log} 2>&1'

checkpoint length_filtration:
    input:
        f'{out_dir}/length_filtration/bmge/{{gene}}.bmge'
    output:
        f'{out_dir}/length_filtration/bmge/{{gene}}.length_filtered'
    params:
        threshold=0.5
    log:
        f'{out_dir}/logs/length_filtration/{{gene}}.log'
    run:
        original_name = f'{wildcards.gene}.length_filtered'
        length = None
        with open(output[0], 'w') as outfile, open(log[0], 'w') as logfile:
            for record in SeqIO.parse(input[0], 'fasta'):
                if length is None:
                    length = len(record.seq)
                coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
                if coverage > params.threshold:
                    outfile.write(f'>{record.description}\n{record.seq}\n')
                else:
                    logfile.write(f'deleted: {record.name} {coverage}')
            
            if os.stat(input[0]).st_size == 0:
                logfile.write(f'All sequences were removed during length filtration')

rule mafft:
    input:
        f'{out_dir}/length_filtration/bmge/{{gene}}.length_filtered'
    output:
        f'{out_dir}/mafft/{{gene}}.aln'
    log:
        f'{out_dir}/logs/mafft/{{gene}}.log'
    shell:
        'mafft --thread 1 --globalpair --maxiterate 1000 --unalignlevel 0.6 {input} >{output} 2>{log}'

rule divvier:
    input:
        f'{out_dir}/mafft/{{gene}}.aln'
    output:
        f'{out_dir}/divvier/{{gene}}.aln.partial.fas',
        f'{out_dir}/divvier/{{gene}}.aln.PP'
    log:
        f'{out_dir}/logs/divvier/{{gene}}.log'
    shell:
        f'''
        divvier -minicol 4 -partial {{input}} >{{log}} 2>{{log}}

        mv {out_dir}/mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/divvier >{{log}} 2>{{log}}
        mv {out_dir}/mafft/{{wildcards.gene}}.aln.PP {out_dir}/divvier >{{log}} 2>{{log}}
        '''

rule trimal:
        input:
            f'{out_dir}/divvier/{{gene}}.aln.partial.fas'
        output:
            f'{out_dir}/trimal/{{gene}}.final'
        log:
            f'{out_dir}/logs/trimal/{{gene}}.log'
        run:
            shell('trimal -in {input} -gt 0.01 -out {output} >{log} 2>{log}')

            records = []
            for record in SeqIO.parse(output[0], 'fasta'):
                if len(str(record.seq).replace('-', '').replace('X','')) > 0:
                    records.append(record)

            SeqIO.write(records, output[0], "fasta")

def get_raxml_input(wildcards):
    gene = '{wildcards.gene}'.format(wildcards=wildcards)
    if trees_only:
        return f'{in_dir}/{gene}.fas'
    else:
        return f'{out_dir}/trimal/{gene}.final'

rule raxml:
    input:
        get_raxml_input
    output:
        f'{out_dir}/raxml/RAxML_bipartitions.{{gene}}.tre',
        f'{out_dir}/raxml/RAxML_bootstrap.{{gene}}.tre',
        f'{out_dir}/raxml/RAxML_bestTree.{{gene}}.tre',
        f'{out_dir}/raxml/RAxML_info.{{gene}}.tre',
        f'{out_dir}/raxml/RAxML_bipartitionsBranchLabels.{{gene}}.tre',
    log:
        f'{out_dir}/logs/raxml/{{gene}}.log'
    params:
        raxml_out=f'{out_dir}/raxml'
    shell:
        'raxmlHPC-PTHREADS-AVX2 -T 1 -m PROTGAMMALG4XF -f a -s {input} -n {wildcards.gene}.tre -w {params.raxml_out} -x 123 -N 100 -p 12345 >{log} 2>{log}'
 