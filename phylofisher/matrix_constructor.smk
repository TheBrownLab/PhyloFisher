import csv
from glob import glob
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

out_dir = config['out_dir']
in_dir = config['in_dir']
in_format = config['in_format']
out_format = config['out_format']
genes = config['genes'].split(',')
concatenation_only = config['concatenation_only']
trimal_gt = config['trimal_gt']

# Accepted out formats with respective suffix
out_dict = {'fasta':          'fas',
            'phylip':         'phy',
            'phylip-relaxed': 'phy',
            'nexus':          'nex'}

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
        'prequal {input} >{log} 2>{log}'


rule mafft:
    input:
        f'{out_dir}/prequal/{{gene}}.aa.filtered'
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
        divvier -minicol 4 -partial -divvygap {{input}} >{{log}} 2>{{log}}

        mv {out_dir}/mafft/{{wildcards.gene}}.aln.partial.fas {out_dir}/divvier >{{log}} 2>{{log}}
        mv {out_dir}/mafft/{{wildcards.gene}}.aln.PP {out_dir}/divvier >{{log}} 2>{{log}}
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
        params:
            gt=trimal_gt
        shell:
            'trimal -in {input} -gt {params.gt} -out {output} >{log} 2>{log}'


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


def get_orgs(input_folder):
    """

    :param input_folder:
    :return:
    """
    name_set = set()
    files = sorted(glob(f'{input_folder}/*'))
    for file in files:
        with open(file) as f:
            for record in SeqIO.parse(f, in_format):
                fname = record.description
                name = fname.split('_')[0]
                name_set.add(name)
    return sorted(list(name_set))


rule construct_matrix:
    input:
        expand(f'{out_dir}/trimal/{{gene}}.final', gene=genes)
    output:
        f'{out_dir}/indices.tsv',
        f'{out_dir}/matrix.{out_dict[out_format.lower()]}',
        f'{out_dir}/matrix_constructor_stats.tsv'
    run:
        with open(output[0], 'w') as outfile:
            outfile.write('Gene\tStart\tStop\n')
            files = sorted(glob(f'{out_dir}/trimal/*.final'))

            total_len = 0
            res_dict = defaultdict(str)
            for file in files:
                gene = os.path.basename(file).split('.')[0]
                length = 0
                seq_dict = {}
                if concatenation_only:
                    myformat = in_format
                else:
                    myformat = 'fasta'
                for record in SeqIO.parse(file, myformat):
                    length = len(record.seq)
                    seq_dict[record.id.split('_')[0]] = str(record.seq)
                start_len = total_len + 1
                total_len += length
                outfile.write(f'{gene}\t{start_len}\t{total_len}\n')
                for org in get_orgs(in_dir):
                    if org in seq_dict:
                        res_dict[org] += seq_dict[org]
                    else:
                        res_dict[org] += ('-' * length)
        
        records = []
        for org, seq in res_dict.items():
            records.append(SeqRecord(Seq(seq),
                                    id=org,
                                    name='',
                                    description=''))
                                    
        if out_format.lower() in out_dict:
            with open(output[1], "w") as handle:
                SeqIO.write(records, handle, out_format.lower())
        else:
            sys.exit('Invalid Output Format')

        with open(output[2], 'w') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow(['Taxon', 'PercentMissingData'])
                missing = []
                for record in SeqIO.parse(input[0], out_format.lower()):
                    missing.append((record.name, (record.seq.count('-') / total_len) * 100))
                for org_missing in sorted(missing, key=lambda x: x[1], reverse=True):
                    tsv_writer.writerow(list(org_missing))

            