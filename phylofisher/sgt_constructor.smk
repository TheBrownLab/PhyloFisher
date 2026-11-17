# The main entry point of your workflow.
import os

from Bio import SeqIO

out_dir = config['out_dir']
in_dir = config['in_dir']
genes = config['genes'].split(',')
trees_only = config['trees_only']
no_trees = config['no_trees']
length_filter_only = config['length_filter_only']
tree_colors = config['tree_colors']
metadata = config['metadata']
input_metadata = config['input_metadata']
threads_per_job = config['threads_per_job']


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
        'prequal {input} &>{log}'

rule length_filter_mafft:
    input:
        f'{out_dir}/prequal/{{gene}}.aa.filtered'
    output:
        f'{out_dir}/length_filter_mafft/{{gene}}.aln'
    log:
        f'{out_dir}/logs/length_filter_mafft/{{gene}}.log'
    conda:
        'mafft.yaml'
    threads:
        threads_per_job
    shell:
        '''
        mafft --thread {threads} --globalpair --maxiterate 1000 --unalignlevel 0.6 {input} >{output} 2>{log}
        [ -s {output} ] || {{ echo "MAFFT produced empty output" >>{log}; exit 1; }}
        '''

rule length_filter_divvier:
    input:
        f'{out_dir}/length_filter_mafft/{{gene}}.aln'
    output:
        f'{out_dir}/length_filter_divvier/{{gene}}.aln.partial.fas',
        f'{out_dir}/length_filter_divvier/{{gene}}.aln.PP'
    log:
        f'{out_dir}/logs/length_filter_divvier/{{gene}}.log'
    conda:
        'divvier.yaml'
    shell:
        f'''
        divvier -mincol 4 -partial {{input}} &>{{log}}

        mv {out_dir}/length_filter_mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/length_filter_divvier >>{{log}} 2>&1
        mv {out_dir}/length_filter_mafft/{{wildcards.gene}}.aln.PP {out_dir}/length_filter_divvier >>{{log}} 2>&1
        '''

rule x_to_dash:
    input:
        f'{out_dir}/length_filter_divvier/{{gene}}.aln.partial.fas'
    output:
        f'{out_dir}/length_filter_bmge/{{gene}}.pre_bmge'
    log:
        f'{out_dir}/logs/x_to_dash/{{gene}}.log'
    run:
        with open(output[0], 'w') as res:
            for record in SeqIO.parse(input[0], 'fasta'):
                res.write(f'>{record.name}\n{str(record.seq).replace("X", "-")}\n')

rule length_filter_bmge:
    input:
        f'{out_dir}/length_filter_bmge/{{gene}}.pre_bmge'
    output:
        f'{out_dir}/length_filter_bmge/{{gene}}.bmge'
    log:
        f'{out_dir}/logs/length_filter_bmge/{{gene}}.log'
    conda:
        'bmge.yaml'
    shell:
        'bmge -t AA -g 0.3 -i {input} -of {output} &>{log}'

checkpoint length_filtration:
    input:
        f'{out_dir}/prequal/{{gene}}.aa',
        f'{out_dir}/length_filter_bmge/{{gene}}.bmge'
    output:
        f'{out_dir}/length_filtered/{{gene}}.length_filtered'
    params:
        threshold=0.5
    log:
        f'{out_dir}/logs/length_filtered/{{gene}}.log'
    run:
        with open(log[0], 'w') as logfile:
            ids_to_keep = set()
            for record in SeqIO.parse(input[1], 'fasta'):
                coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
                if coverage > params.threshold:
                    ids_to_keep.add(record.description)
                else:
                    logfile.write(f'deleted: {record.name} {coverage}\n')

            if not ids_to_keep:
                logfile.write('All sequences were removed during length filtration\n')

            # Write the unaligned sequences that passed the coverage filter
            with open(output[0], 'w') as outfile:
                for record in SeqIO.parse(input[0], 'fasta'):
                    if record.description in ids_to_keep:
                        outfile.write(f'>{record.description}\n{str(record.seq)}\n')


rule mafft:
    input:
        f'{out_dir}/length_filtered/{{gene}}.length_filtered'
    output:
        f'{out_dir}/mafft/{{gene}}.aln'
    log:
        f'{out_dir}/logs/mafft/{{gene}}.log'
    conda:
        'mafft.yaml'
    threads:
        threads_per_job
    shell:
        '''
        mafft --thread {threads} --globalpair --maxiterate 1000 --unalignlevel 0.6 {input} >{output} 2>{log}
        [ -s {output} ] || {{ echo "MAFFT produced empty output" >>{log}; exit 1; }}
        '''

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
        divvier -mincol 4 -partial {{input}} &>{{log}}

        mv {out_dir}/mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/divvier >>{{log}} 2>&1
        mv {out_dir}/mafft/{{wildcards.gene}}.aln.PP {out_dir}/divvier >>{{log}} 2>&1
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
        'trimal -in {input} -gt 0.01 -out {output} &>{log}'


rule remove_gaps:
    input:
        f'{out_dir}/trimal/{{gene}}.trimal'
    output:
        f'{out_dir}/trimal/{{gene}}.final'
    run:
        records = []
        for record in SeqIO.parse(input[0], 'fasta'):
            if len(str(record.seq).replace('-', '').replace('X', '')) > 0:
                records.append(record)

        SeqIO.write(records, output[0], 'fasta')


def get_gene_done_input(wildcards):
    '''
    Per-gene outputs required before a gene is considered done.
    Waits only on this gene's length_filtration checkpoint, so each gene moves
    through the pipeline independently. Genes emptied by length filtration
    require nothing and are skipped.
    '''
    gene = wildcards.gene
    if trees_only:
        return [f'{out_dir}-local/trees/{gene}.final',
                f'{out_dir}-local/trees/{gene}.treefile']

    length_filtered = checkpoints.length_filtration.get(gene=gene).output[0]
    if os.path.getsize(length_filtered) == 0:
        return []

    if no_trees:
        return [f'{out_dir}/trimal/{gene}.final']

    return [f'{out_dir}-local/trees/{gene}.trimmed',
            f'{out_dir}-local/trees/{gene}.untrimmed',
            f'{out_dir}-local/trees/{gene}.final',
            f'{out_dir}-local/trees/{gene}.treefile']

rule gene_done:
    input:
        get_gene_done_input
    output:
        touch(f'{out_dir}/done/{{gene}}.done')

rule no_trees_done:
    input:
        [f'{out_dir}/done/{gene}.done' for gene in genes]
    output:
        touch(f'{out_dir}/trimal/no_trees.done')


def get_raxml_input(wildcards):
    if trees_only:
        return f'{in_dir}/{wildcards.gene}.fas'
    else:
        return f'{out_dir}/trimal/{wildcards.gene}.final'

if not no_trees:
    rule raxml:
        input:
            get_raxml_input
        output:
            f'{out_dir}/raxml/{{gene}}.raxml.rba',
            f'{out_dir}/raxml/{{gene}}.raxml.startTree',
            f'{out_dir}/raxml/{{gene}}.raxml.bestTree',
            f'{out_dir}/raxml/{{gene}}.raxml.mlTrees',
            f'{out_dir}/raxml/{{gene}}.raxml.support',
            f'{out_dir}/raxml/{{gene}}.raxml.bestModel',
            f'{out_dir}/raxml/{{gene}}.raxml.bootstraps',
            f'{out_dir}/raxml/{{gene}}.raxml.log'
        log:
            f'{out_dir}/logs/raxml/{{gene}}.log'
        conda:
            'raxml-ng.yaml'
        params:
            raxml_out=f'{out_dir}/raxml'
        threads:
            threads_per_job
        shell:
            '''
            raxml-ng --all \
                --msa {input} \
                --prefix {params.raxml_out}/{wildcards.gene} \
                --model LG4X+G4 \
                --tree pars{{10}} \
                --bs-trees 100 \
                --threads {threads} &>{log}
            '''

if not no_trees and trees_only:
    rule cp_trees:
        input:
            f'{in_dir}/{{gene}}.fas',
            f'{out_dir}/raxml/{{gene}}.raxml.support'
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

elif not no_trees:
    rule cp_trees:
        input:
            f'{out_dir}/length_filter_bmge/{{gene}}.bmge',
            f'{out_dir}/trimal/{{gene}}.final',
            f'{out_dir}/raxml/{{gene}}.raxml.support',
            f'{out_dir}/mafft/{{gene}}.aln'
        output:
            f'{out_dir}/trees/{{gene}}.trimmed',
            f'{out_dir}/trees/{{gene}}.final',
            f'{out_dir}/trees/{{gene}}.treefile',
            f'{out_dir}-local/trees/{{gene}}.trimmed',
            f'{out_dir}-local/trees/{{gene}}.final',
            f'{out_dir}-local/trees/{{gene}}.treefile',
            f'{out_dir}-local/trees/{{gene}}.untrimmed'
        shell:
            '''
            cp {input[0]} {output[0]}
            cp {input[0]} {output[3]}
            cp {input[1]} {output[1]}
            cp {input[1]} {output[4]}
            cp {input[2]} {output[2]}
            cp {input[2]} {output[5]}
            cp {input[3]} {output[6]}
            '''

rule cp_metadata:
    input:
        f'{tree_colors}',
        f'{metadata}',
        f'{input_metadata}',
    output:
        f'{out_dir}-local/tree_colors.tsv',
        f'{out_dir}-local/metadata.tsv',
        f'{out_dir}-local/input_metadata.tsv'
    shell:
        '''
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        cp {input[2]} {output[2]}
        '''

def get_tar_local_dir_input(wildcards):
    ret = []

    if length_filter_only:
        # Only package length-filtered sequences
        for gene in genes:
            ret.append(f'{out_dir}/length_filtered/{gene}.length_filtered')
    else:
        # Per-gene outputs are pulled in by gene_done; genes emptied by length filtration are skipped
        for gene in genes:
            ret.append(f'{out_dir}/done/{gene}.done')

        ret.append(f'{out_dir}-local/tree_colors.tsv')
        ret.append(f'{out_dir}-local/metadata.tsv')
        ret.append(f'{out_dir}-local/input_metadata.tsv')

    return ret

if not length_filter_only:
    rule tar_local_dir:
        input:
            get_tar_local_dir_input
        output:
            f'{out_dir}-local.tar.gz'
        log:
            f'{out_dir}/logs/tar_local_dir.log'
        params:
            out_dir_parent=os.path.dirname(out_dir),
            out_dir_base=f'{os.path.basename(out_dir)}-local'
        shell:
            f'''
            tar -czvf {out_dir}-local.tar.gz -C {{params.out_dir_parent}} {{params.out_dir_base}} &>{{log}}
            rm -r {out_dir}-local
            '''
else:
    # When length_filter_only is set, create a tarball of length-filtered sequences
    rule tar_length_filtered:
        input:
            get_tar_local_dir_input
        output:
            f'{out_dir}-length_filtered.tar.gz'
        log:
            f'{out_dir}/logs/tar_length_filtered.log'
        shell:
            f'''
            tar -czvf {out_dir}-length_filtered.tar.gz -C {out_dir} length_filtered &>{{log}}
            '''
