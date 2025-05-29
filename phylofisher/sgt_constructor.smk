# The main entry point of your workflow.
import os
from peewee import *
from phylofisher.db_map import database, BaseModel, Genes, Taxonomies, Metadata, Sequences
from Bio import SeqIO

out_dir = config['out_dir']
in_dir = config['in_dir']
genes = config['genes'].split(',')
trees_only = config['trees_only']
no_trees = config['no_trees']
pf_database = config['database']
input_metadata = config['input_metadata']


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
    conda:
        'prequal.yaml'
    shell:
        'prequal {input} &> {log}'

rule length_filter_mafft:
    input:
        f'{out_dir}/prequal/{{gene}}.aa.filtered'
    output:
        f'{out_dir}/length_filtration/mafft/{{gene}}.aln'
    log:
        f'{out_dir}/logs/length_filter_mafft/{{gene}}.log'
    conda:
        'mafft.yaml'
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
    conda:
        'divvier.yaml'
    shell:
        f'''
        divvier -mincol 4 -partial {{input}} >{{log}} 2>{{log}}

        mv {out_dir}/length_filtration/mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/length_filtration/divvier &> {{log}}
        mv {out_dir}/length_filtration/mafft/{{wildcards.gene}}.aln.PP {out_dir}/length_filtration/divvier &> {{log}}
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
    conda:
        'bmge.yaml'
    shell:
        'bmge -t AA -g 0.3 -i {input} -of {output} &> {log}'

rule length_filtration:
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
    conda:
        'mafft.yaml'
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
    conda:
        'divvier.yaml'
    shell:
        f'''
        divvier -minicol 4 -partial {{input}} >{{log}} 2>{{log}}

        mv {out_dir}/mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/divvier &>{{log}}
        mv {out_dir}/mafft/{{wildcards.gene}}.aln.PP {out_dir}/divvier &>{{log}}
        '''

rule trimal:
        input:
            f'{out_dir}/divvier/{{gene}}.aln.partial.fas'
        output:
            f'{out_dir}/trimal/{{gene}}.trimal'
        log:
            f'{out_dir}/logs/trimal/{{gene}}.log'
        conda:
            'trimal.yaml'
        shell:
            'trimal -in {input} -gt 0.01 -out {output} &> {log}'


rule remove_gaps:
        input:
            f'{out_dir}/trimal/{{gene}}.trimal'
        output:
            f'{out_dir}/trimal/{{gene}}.final'
        log:
            f'{out_dir}/logs/trimal/{{gene}}.log'
        run:
            records = []
            for record in SeqIO.parse(input[0], 'fasta'):
                if len(str(record.seq).replace('-', '').replace('X','')) > 0:
                    records.append(record)

            SeqIO.write(records, output[0], "fasta")
            

def get_raxml_input(wildcards):
    gene = '{wildcards.gene}'.format(wildcards=wildcards)
    if trees_only:
        return f'{in_dir}/{gene}.fas'
    else:
        return f'{out_dir}/trimal/{gene}.final'

rule iqtree:
    input:
        get_raxml_input
    output:
        f'{out_dir}/iqtree/{{gene}}.treefile'
    log:
        f'{out_dir}/logs/iqtree/{{gene}}.log'
    conda:
        'iqtree.yaml'
    params:
        iqtree_out=lambda wildcards: f'{out_dir}/iqtree/{wildcards.gene}'
    shell:
        '''
        iqtree -s {input} \
        -pre {params.iqtree_out} \
        -m ELM+C20 \
        -B 1000 \
        -alrt 1000 \
        -T 1 \
        &> {log}
        '''

if trees_only:
    rule cp_trees:
        input:
            f'{in_dir}/{{gene}}.fas',
            f'{out_dir}/iqtree/{{gene}}.treefile'
        output:
            f'{out_dir}/trees/{{gene}}.final',
            f'{out_dir}/trees/{{gene}}.treefile',
            f'{out_dir}-local/trees/{{gene}}.final',
            f'{out_dir}-local/trees/{{gene}}.treefile'
        shell:
            '''
            cp {input[0]} {output[0]}
            cp {input[0]} {output[2]}
            cp {input[1]} {output[1]}
            cp {input[1]} {output[3]}
            '''
            
else:
    rule cp_trees:
        input:
            f'{out_dir}/length_filtration/bmge/{{gene}}.length_filtered',
            f'{out_dir}/trimal/{{gene}}.final',
            f'{out_dir}/iqtree/{{gene}}.treefile'
        output:
            f'{out_dir}/trees/{{gene}}.trimmed',
            f'{out_dir}/trees/{{gene}}.final',
            f'{out_dir}/trees/{{gene}}.treefile',
            f'{out_dir}-local/trees/{{gene}}.trimmed',
            f'{out_dir}-local/trees/{{gene}}.final',
            f'{out_dir}-local/trees/{{gene}}.treefile'
        shell:
            '''
            cp {input[0]} {output[0]}
            cp {input[0]} {output[3]}
            cp {input[1]} {output[1]}
            cp {input[1]} {output[4]}
            cp {input[2]} {output[2]}
            cp {input[2]} {output[5]}
            '''

rule cp_input_metadata:
    input:
        f'{input_metadata}',
    output:
        f'{out_dir}-local/input_metadata.tsv'
    shell:
        '''
        cp {input[0]} {output[0]}
        '''

rule make_metadata_tsv:
    input:
        f'{pf_database}',
    output:
        f'{out_dir}-local/metadata.tsv',
        f'{out_dir}-local/tree_colors.tsv'
    run:
        database.init(input[0])
        database.connect()
        db_query = Metadata.select()
        tax_query = Taxonomies.select()

        with open(output[0], 'w') as meta_out, open(output[1], 'w') as color_out:
            meta_out.write('Unique ID\tLong Name\tHigher Taxonomy\tLower Taxonomy\tData Type\tSource\n')
            color_out.write('Taxonomy\tColor\n')

            tax_dict = {}
            for taxa in tax_query:
                if taxa.color:
                    color_out.write(f'{taxa.taxonomy}\t{taxa.color}\n')
                tax_dict[taxa.id] = taxa.taxonomy

            for q in db_query:
                meta_out.write(f'{q.short_name}\t{q.long_name}\t{tax_dict[int(str(q.higher_taxonomy))]}\t{tax_dict[int(str(q.lower_taxonomy))]}\t{q.data_type}\t{q.source}\n')
            
        database.close()

def get_tar_local_dir_input(wildcards):
    ret = []
    for gene in genes:
        if not trees_only:
            ret.append(f'{out_dir}-local/trees/{gene}.trimmed')
        ret.append(f'{out_dir}-local/trees/{gene}.final')
        ret.append(f'{out_dir}-local/trees/{gene}.treefile')

    ret.append(f'{out_dir}-local/tree_colors.tsv')
    ret.append(f'{out_dir}-local/metadata.tsv')
    ret.append(f'{out_dir}-local/input_metadata.tsv')

    return ret

rule tar_local_dir:
    input:
        get_tar_local_dir_input
    output:
        f'{out_dir}-local.tar.gz'
    log:
        f'{out_dir}/logs/tar_local_dir.log'
    params:
        out_dir_base=f'{out_dir.split("/")[-1]}-local'
    shell:
        f'''
        tar -czvf {out_dir}-local.tar.gz {{params.out_dir_base}} >{{log}} 2>{{log}}
        rm -r {out_dir}-local
        '''
