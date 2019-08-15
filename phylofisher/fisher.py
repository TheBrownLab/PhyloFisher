#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import subprocess
import argparse
from multiprocessing import Pool
from functools import partial
from pathlib import Path
import configparser
import warnings
from shutil import copyfile, rmtree
from Bio import SeqIO
from ete3 import Tree
from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO


class Hit:

    def __init__(self, name, seq, query, hmm_hits):
        self.name = name
        self.seq = seq
        self.query = query
        self.hmm_hits = hmm_hits
        self.path_ = ""
        self.quality = ""

    def in_hmm_hits(self):
        if self.name in self.hmm_hits:
            return True
        return False


class Query:

    def __init__(self, query, hmm_hits, infile_proteins):
        self.query = query
        self.hmm_hits = hmm_hits
        self.infile_proteins = infile_proteins
        self.squery = False

    def hits(self):
        for hit in self.hmm_hits:
            hit = Hit(hit, self.infile_proteins[hit], self.query, self.hmm_hits)
            hit.path_ += "HMM"
            yield hit


class SpecQuery:

    def __init__(self, query, spec_queries, hmm_hits, infile_proteins):
        self.query = query
        self.hmm_hits = hmm_hits
        self.seq = None
        self.infile_proteins = infile_proteins
        self.organisms = [org for org in spec_queries.split(',')]

    def get_specific_query(self):
        gene_dict = {}
        for record in SeqIO.parse(str(Path(dfo, f'orthologs/{self.query}.fas')), 'fasta'):
            gene_dict[record.name] = str(record.seq)
        for org in self.organisms:
            if org in gene_dict:
                self.seq = gene_dict[org]
                return org

    def spec_query_blast(self):
        db = f"tmp/{sample_name}/{os.path.basename(infile)}.blastdb"
        qfile = f'tmp/{sample_name}/{self.query}.fas'
        bout = f'tmp/{sample_name}/{self.query}.blastout'
        with open(qfile, 'w') as f:
            f.write(f'>{self.query}\n{self.seq.replace("-", "")}')
        cmd = f'blastp -evalue 1e-10 -query {qfile} -db {db} -out {bout} -outfmt "6 qseqid sseqid evalue"'
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        blast_hits = []
        for line_ in open(bout):
            prot_name = line_.split('\t')[1]
            if prot_name not in blast_hits:
                blast_hits.append(prot_name)
        return blast_hits

    def hits(self):
        spec_query = self.get_specific_query()
        if spec_query:
            blast_hits = self.spec_query_blast()
            if not blast_hits:
                return None
            for blast_hit in blast_hits:
                hit = Hit(blast_hit, self.infile_proteins[blast_hit], self.query, self.hmm_hits)
                hit.path_ += "SBH"
                yield hit
        else:
            for hmm_hit in self.hmm_hits:
                hit = Hit(hmm_hit, self.infile_proteins[hmm_hit], self.query, self.hmm_hits)
                hit.path_ += "HMM"
                yield hit


def length_check(trimmed_aln):
    correct_length = {}
    for record in SeqIO.parse(trimmed_aln, 'fasta'):
        if '@' in record.name:
            true_length = len(str(record.seq).replace('-', '').replace('X', ''))
            if true_length / len(record.seq) > 0.3:
                correct_length[record.name] = round((true_length / len(record.seq)), 2)
    return correct_length


def makeblastdb():
    cmd = f"makeblastdb -in {infile} -out tmp/{sample_name}/{os.path.basename(infile)}.blastdb -dbtype prot"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)


def hmmer(query):
    hmm_prof = str(Path(dfo, f'profiles/{query}.hmm'))
    cmd = f"hmmsearch -E 1e-10 {hmm_prof} {infile} > tmp/{sample_name}/{query}.hmmout"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
    hits = []
    for hit in SearchIO.read(f"tmp/{sample_name}/{query}.hmmout", "hmmer3-text"):
        hits.append(hit.id)
    return query, hits


def get_gene_dict(threads, infile_proteins, spec_queries=None):
    gene_dict = {}
    with Pool(processes=threads) as pool:
        hmm_pool = list(pool.map(hmmer, profiles))
        for query, hmm_hits in hmm_pool:
            if hmm_hits:
                if spec_queries:
                    gene_dict[query] = SpecQuery(
                        query, spec_queries, hmm_hits, infile_proteins)
                else:
                    gene_dict[query] = Query(query, hmm_hits, infile_proteins)
            else:
                print(f"{query} exluded. No hmm hits.")
    return gene_dict


def best_hits(max_hits, gene):
    candidates = []
    hits = gene.hits()
    counter = 0
    if hits:
        for hit in hits:
            if counter == max_hits:
                break
            counter += 1
            if hit.in_hmm_hits() is True:
                candidates.append(hit)
    return gene.query, candidates


def bac_gog_db():
    bac = set()
    for line_ in open(str(Path(dfo, 'orthomcl/bacterial'))):
        bac.add(line_[:-1])
    g_og = {}
    for line_ in open(str(Path(dfo, 'orthomcl/gene_og'))):
        sline = line_.split('\t')
        g_og[sline[0]] = sline[1]
    return bac, g_og


def get_infile_proteins():
    infile_proteins = {}
    for record in SeqIO.parse(infile, 'fasta'):
        infile_proteins[record.id] = str(record.seq)
    return infile_proteins


def get_candidates(threads, max_hits, queries):
    with Pool(processes=threads) as pool:
        func = partial(best_hits, max_hits)
        candidates = list(pool.map(func, queries))
        return candidates


def get_hmm_profiles():
    hmm_profiles = []
    for file in os.listdir(str(Path(dfo, 'profiles/'))):
        hmm_profiles.append(file.split('.')[0])
    return hmm_profiles


def makedirs():
    directories = ['tmp', 'fasta']
    for directory in directories:
        if not os.path.isdir(directory):
            os.makedirs(directory)


def phylofisher(threads, max_hits, spec_queries=None):
    infile_proteins = get_infile_proteins()
    if spec_queries:
        makeblastdb()
        gene_dict = get_gene_dict(threads, infile_proteins, spec_queries)
    else:
        gene_dict = get_gene_dict(threads, infile_proteins)

    candidates = get_candidates(threads, max_hits, gene_dict.values())

    with open('tmp/for_diamond.fasta', 'a') as f:
        for gene, candi_list in candidates:
            if candi_list:
                for hit in candi_list:
                    f.write(f'>{hit.name}_{hit.path_}@{gene}\n{hit.seq}\n')


def taxonomy_dict():
    tax_g = {}
    for line_ in open(str(Path(dfo, 'metadata.tsv'))):
        if 'Full Name' not in line_:
            sline = line_.split('\t')
            tax = sline[0].strip()
            group = sline[2].strip()
            tax_g[tax] = group
    return tax_g


def cluster_rename_sequences():
    clustered = f'tmp/{sample_name}/clustered.fasta'
    cmd = f'cd-hit -i {fasta_file} -o {clustered} -c 0.98'
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    original_names = {}
    n = 1
    with open(f'tmp/{sample_name}/clustered_renamed.fasta', 'w') as res:
        for record in SeqIO.parse(clustered, 'fasta'):
            new_name = f"{sample_name}_{n}"
            original_names[new_name] = record.name
            res.write(f'>{new_name}\n{record.seq}\n')
            n += 1
    abs_path = os.path.abspath(f'tmp/{sample_name}/clustered_renamed.fasta')
    with open('original_names.tsv', 'a') as f:
        for new_name, or_name in original_names.items():
            f.write(f'{or_name}\t{new_name}\n')
    return abs_path


def check_input():
    #TODO check is file is aa not nucl
    errors = ''
    taxonomic_groups = tax_group.values()
    n = 1
    for line in open(multi_input):
        if "FILE_NAME" not in line:
            n += 1
            metadata_input = line.split('\t')
            f_file = str(Path(metadata_input[0], metadata_input[1]))
            if os.path.isfile(f_file.strip()) is False:
                errors += f"line {n}: file {f_file} doesn't exist\n"
            s_name = metadata_input[2].strip()
            if s_name  in tax_group:
                errors += f'line {n}: {s_name} already in metadata\n'
            s_queries = metadata_input[5]
            if s_queries.lower().strip() != 'none':
                for q in s_queries.split(','):
                    if q.strip() not in tax_group:
                        errors += f'line {n}: {q} not in metadata\n'
                tax = metadata_input[3].strip()
                if '*' not in tax:
                    if tax not in taxonomic_groups:
                        errors += f'line {n}: {tax} not in metadata\n'
    if errors:
        sys.exit(f'Please check your input file:\n{errors[:-1]}')



def diamond():
    for_diamond = 'tmp/for_diamond.fasta'

    db = str(Path(dfo, 'orthomcl/orthomcl.diamonddb'))
    out = 'tmp/diamond.res'
    cmd = (
        f'diamond blastp -e 1e-10 -q {for_diamond} --more-sensitive '
        f'--db {db} -o {out} -p {args.threads} --outfmt 6 qseqid stitle evalue')
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    datasetdb = str(Path(dfo, 'datasetdb/datasetdb.dmnd'))
    out2 = 'tmp/dataset_diamond.res'
    cmd2 = (
        f'diamond blastp -e 1e-10 -q {for_diamond} --more-sensitive '
        f'--db {datasetdb} -o {out2} -p {args.threads} --outfmt 6 qseqid stitle evalue')
    subprocess.run(cmd2, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def correct_phylo_group(parent, sample_taxomomy):
    # TODO add paralog filter
    groups = set()
    for org in parent.get_leaf_names():
        if org in tax_group:
            groups.add(tax_group[org])
    if len(groups) > 2:
        return False
    elif sample_taxomomy in groups:
        return True
    return correct_phylo_group(parent.up, sample_taxomomy)


def fasttree(checked_hits):
    full_name = checked_hits[0].name
    org = full_name.split('_')[0]
    gene = full_name.split('@')[1]
    fas = f'tmp/{org}/{gene}.for_ftree'
    aln = f'tmp/{org}/{gene}.aln'
    trim = f'tmp/{org}/{gene}.bmge'
    tree_file = f'tmp/{org}/{gene}.tree'
    copyfile(str(Path(dfo, f'orthologs/{gene}.fas')), f'tmp/{org}/{gene}.for_ftree')
    with open(fas, 'a') as f:
        for hit in checked_hits:
            f.write(f'>{hit.name}\n{hit.seq}\n')
    cmd1 = f'mafft --auto --reorder {fas} > {aln}'
    subprocess.run(cmd1, shell=True, stderr=subprocess.DEVNULL)
    cmd2 = f"java -jar {bmge} -t AA -m BLOSUM30 -b 2 -g 0.6 -i {aln} -of {trim}"
    subprocess.run(cmd2, shell=True, stdout=subprocess.DEVNULL)
    cmd3 = f"fasttree {trim} > {tree_file}"
    subprocess.run(cmd3, shell=True, stderr=subprocess.DEVNULL)
    tree = Tree(tree_file)
    correct_len = length_check(trim)
    good_hits = []
    bb_hits = []
    for hit in checked_hits:
        if hit.name in correct_len:
            bb_hits.append(hit)
            hit_node = tree.search_nodes(name=hit.name)[0]
            if correct_phylo_group(hit_node.up, input_taxonomy[org]) is True:
                good_hits.append(hit)
    if len(good_hits) > 0:
        return good_hits
    elif len(bb_hits) > 0:
        for hit in bb_hits:
            hit.name = hit.name.replace('_SBH', '_BBH')
        return bb_hits
    else:
        return []


def parse_diamond_output():
    correct_hits = set()
    proccesed = set()
    for line in open('tmp/diamond.res'):
        sline = line.split("\t")
        full_name = sline[0]
        hit, gene = full_name.split('@')
        if full_name not in proccesed:
            proccesed.add(full_name)
            org = sline[1].split('|')[0]
            og = sline[1].split('|')[2].strip()
            #TODO are we sure that we want to control only best hit?
            if org not in bacterial and og in gene_og[gene]:
                correct_hits.add(full_name)
    return correct_hits


def new_best_hits(candidate_hits):
    top_candidates = []
    if candidate_hits:
        if '_SBH' in candidate_hits[0].name:
            top_candidates = fasttree(candidate_hits)
        else:
            top_candidates = candidate_hits

    if top_candidates:
        gene = top_candidates[0].name.split('@')[1]
        dataset = f'fasta/{gene}.fas'
        if not os.path.isfile(dataset):
            copyfile(str(Path(dfo, f'orthologs/{gene}.fas')), dataset)
        n = 0
        for cand in top_candidates:
            n += 1
            with open(dataset, 'a') as d:
                seq_name = cand.name.split("@")[0]
                gene = cand.name.split("@")[1]
                try:
                    if gene == reciprocal_hits[seq_name[:-4]]:
                        d.write(f'>{seq_name}_q{n}r\n{cand.seq}\n')
                    else:
                        d.write(f'>{seq_name}_q{n}n\n{cand.seq}\n')
                        with open('nonreciprocal_hits.txt', 'a') as nonrep:
                            nonrep.write(f'nonreciprocal hit:{cand.name}; Best hit from:{reciprocal_hits[seq_name[:-4]]}\n')
                            print(f'nonreciprocal hit:{cand.name}; Best hit from: {reciprocal_hits[seq_name[:-4]]}')
                except KeyError:
                    n -= 1


def main_func(gene_hits):
    with Pool(processes=int(args.threads)) as pool:
        pool.map(new_best_hits, gene_hits)


def prepare_good_hits():
    gene_hits = defaultdict(list)
    for record in SeqIO.parse('tmp/for_diamond.fasta', 'fasta'):
        if record.name in correct_hits:
            gene = record.name.split('@')[1]
            org = record.name.split('_')[0]
            org_gene = f'{org}_{gene}'
            gene_hits[org_gene].append(record)
    main_func(list(gene_hits.values()))


def get_reciprocal_hits():
    reciprocal = {}
    proccesed = set()
    for line in open('tmp/dataset_diamond.res'):
        sline = line.split("\t")
        full_name = sline[0]
        sequence, gene = full_name.split('@')
        sequence = sequence[:-4]
        if sequence not in proccesed:
            proccesed.add(sequence)
            hit_gene = sline[1].split('@')[0]
            reciprocal[sequence] = hit_gene
    return reciprocal


def additions_to_input():
    original_input = config['PATHS']['input_file']
    with open(original_input, 'a') as ori:
        for line in open(args.add):
            if "FILE_NAME" not in line:
                ori.write(line)


if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('config.ini')
    parser = argparse.ArgumentParser(description='Script for ortholog fishing.', usage="fisher.py [OPTIONS]")
    parser.add_argument('-t', '--threads', type=int,
                        help='Number of threads, default:1', default=1)
    parser.add_argument('-n', '--max_hits', type=int,
                        help='Max number of hits to check. Default = 1.', default=5)
    parser.add_argument('-v', '--version', action='version', version='0.1')
    parser.add_argument('--keep_tmp', action='store_true')
    parser.add_argument('--add', help='Input file (different from original one in config.ini'
                                      ' only with new organisms. ')
    args = parser.parse_args()

    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())

    tax_group = taxonomy_dict()
    if args.add:
        multi_input = os.path.abspath(args.add)
    else:
        multi_input = os.path.abspath(config['PATHS']['input_file'])
    bmge = config['PATHS']['bmge']
    check_input()

    bacterial, gene_og = bac_gog_db()
    profiles = get_hmm_profiles()

    input_taxonomy = {}

    if not args.add:
        makedirs()
    else:
        if os.path.exists('tmp/'):
            try:
                os.remove('tmp/dataset_diamond.res')
                os.remove('tmp/diamond.res')
                os.remove('tmp/for_diamond.fasta')
            except OSError:
                pass
        else:
            os.mkdir('tmp')
    for line in open(multi_input):
        total_profiles = len(profiles)
        if "FILE_NAME" not in line:
            metadata_input = line.split('\t')
            fasta_file = str(Path(metadata_input[0].strip(), metadata_input[1].strip()))
            sample_name = metadata_input[2].strip()
            taxonomy = metadata_input[3].strip()
            input_taxonomy[sample_name] = taxonomy
            specific_queries = metadata_input[5].strip()
            print(f'{sample_name} has started\n--------------------------')
            os.mkdir(f'tmp/{sample_name}')
            infile = cluster_rename_sequences()
            specific_queries = specific_queries.strip()
            if specific_queries.lower() == 'none':
                specific_queries = None
            phylofisher(args.threads,
                        args.max_hits, specific_queries)
    diamond()
    correct_hits = parse_diamond_output()
    reciprocal_hits = get_reciprocal_hits()
    prepare_good_hits()
    if args.add:
        additions_to_input()
    if not args.keep_tmp:
        rmtree("tmp/")
