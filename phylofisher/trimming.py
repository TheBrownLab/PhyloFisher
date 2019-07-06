def make_alignment(dataset):
    name_root = dataset.split('/')[-1].split('.')[0]
    aln = f'alignments/{name_root}.aln'
    if args.prequal:
        preq = str(Path(dfo, 'lib/prequal-master/prequal'))
        subprocess.run(f'{preq} {dataset}', shell=True, stdout=subprocess.DEVNULL)
        cmd = f'mafft --auto --reorder {dataset}.filtered > {aln}'
    else:
        cmd = f'mafft --auto --reorder {dataset} > {aln}'
    subprocess.run(cmd, shell=True, stderr=subprocess.DEVNULL)
    trimmed = f"trimmed/{name_root}.bmge"
    bmge = str(Path(dfo, 'lib/BMGE-1.12/BMGE.jar'))
    cmd2 = f"java -jar {bmge} -t AA -m BLOSUM30 -b 2 -g 0.6 -i {aln} -of {trimmed}"
    subprocess.run(cmd2, shell=True, stdout=subprocess.DEVNULL)
    good_length(trimmed)


def make_alignments(threads):
    files = glob.glob('fasta/*.fas')
    with Pool(processes=threads) as pool:
        pool.map(make_alignment, files)

def good_length(trimmed_aln):
    records = list(SeqIO.parse(trimmed_aln, 'fasta'))
    with open(trimmed_aln, 'w') as res:
        for record in records:
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > 0.3:
                res.write(f'>{record.description}_{round(coverage, 2)}\n{record.seq}\n')