from Bio import SeqIO
import glob

def good_length(trimmed_aln):
    records = list(SeqIO.parse(trimmed_aln, 'fasta'))
    with open(trimmed_aln, 'w') as res:
        for record in records:
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > 0.3:
                res.write(f'>{record.description}_{round(coverage, 2)}\n{record.seq}\n')
            else:
            	print(record.name, trimmed_aln)


files = glob.glob('*.fas')

for file in files:
	good_length(file)
