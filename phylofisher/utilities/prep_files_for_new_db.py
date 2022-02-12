import sys
import os
import glob
from Bio import SeqIO
import shutil
import argparse

parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')
required.add_argument('-t', '--taxa_list', type = str, help = 'List of taxa as unique IDs to include in new database', required=True)
required.add_argument('-d', '--master_db', type = str, help = 'Path to master phylofisher database', required = True)
required.add_argument('-o', '--out_dir', type = str, help = 'Path to location where output directory will be made', required = True)

args=parser.parse_args()

## Open Taxa list for new DB
infile = open(args.taxa_list)
lines = infile.read()
infile.close

lines = lines.split('\n')
lines=[line for line in lines if line.strip() !=""] #remove empty lines

# get masterDB path
masterdbpath = args.master_db

#Make output directories
outdir = args.out_dir + "/database"
orthooutdir = outdir + "/orthologs"
paraoutdir = outdir + "/paralogs"
protoutdir = outdir + "/proteomes"
try:
	os.mkdir(outdir)
	print(outdir + " created")
except OSError as error:
	print(outdir + " already exists")
	pass
try:
	os.mkdir(orthooutdir)
	print(orthooutdir + " created")
except OSError as error:
	print(orthooutdir + " already exists")
	pass
try:
	os.mkdir(paraoutdir)
	print(paraoutdir + " created")
except OSError as error:
	print(paraoutdir + " already exists")
	pass
try:
	os.mkdir(protoutdir)
	print(protoutdir + " created")
except OSError as error:
	print(protoutdir + " already exists")
	pass
	
## Make new ortholog fastas from selected taxa

orthopath = masterdbpath + "/orthologs/"
orthofastas = (glob.glob(orthopath +"/*.fas"))
for i in orthofastas:
	keep = []
	fasta = SeqIO.parse(open(i), 'fasta')
	fname = i.split('/')[-1]
	print(fname)
	with open(orthooutdir + "/" + fname, "w") as outfasta:
		for record in fasta:
			for line in lines:
				if line in record.id:
					keep.append(record)
					SeqIO.write(record, outfasta, "fasta")
				else:
					pass
					
## Make new paralog fastas from selected taxa
parapath = masterdbpath + "/paralogs/"
parafastas = (glob.glob(parapath +"/*.fas"))
for i in parafastas:
	keep = []
	fasta = SeqIO.parse(open(i), 'fasta')
	fname = i.split('/')[-1]
	print(fname)
	with open(paraoutdir + "/" + fname, "w") as outfasta:
		for record in fasta:
			for line in lines:
				if line in record.id:
					keep.append(record)
					SeqIO.write(record, outfasta, "fasta")
				else:
					pass

## Copy proteome files to new directory
protpath = masterdbpath + "/proteomes/"
proteomes = (glob.glob(protpath +"/*.gz"))
for line in lines:
	shutil.copy(protpath + line + ".faa.tar.gz", protoutdir, follow_symlinks=True)
	print(protpath + line + ".faa.tar.gz")

## Make new metadata file
metain = open(masterdbpath + "/metadata.tsv")
mlines = metain.read()
metain.close

mlines = mlines.split('\n')
mlines=[line for line in mlines if line.strip() !=""]
newmd = ['Unique ID\tLong Name\tHigher Taxonomy\tLower Taxonomy\tData Type\tSource']
for mline in mlines:
	uniqid = mline.split('\t')[0]
	for line in lines:
		if line == uniqid:
			newmd.append(mline)
with open(outdir + '/metadata.tsv', 'w') as metaout:
	for line in newmd:
		metaout.write(line + '\n')
	metaout.close