#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

###NOTE this script is currently very inefficient for processing large datasets. I should add some
### multiprocessing to speed tings up

parser = argparse.ArgumentParser(description="Take Orthofinder output file and produce separate nucleotide fasta files for each orthogroup")

parser.add_argument('--orthogroups', required=True, help='Orthogroups.tsv file produced by Orthofinder', action='store')
parser.add_argument('--fasta_dir', required=True, help='Directory containing nucleotide coding sequence files in fasta format', action='store')
parser.add_argument('--out_dir', required=True, help='Output directory where orthogroup fasta files will be written', action='store')

args=parser.parse_args()

#make the output directory if it doesn't already exist

if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

file =  open(args.orthogroups, 'r')
next(file, None) # skip header

# for each row in the file, get the orthogroup ID and a list of the sequence IDs
for row in file:
	ortho_id = row.split('\t')[0]
	seq_ids = row.split(ortho_id)[1].replace('\t',',')[1:].replace(' ', '').rstrip().split(',')
	
	# create empty list of sequences
	ortho_seqs = []
	
	# Loop over each fasta file in the input directory, check if the record description matches an element in the #
	# sequence IDs list, and if so keep it #
	
	print(
	'''
	==================================
	Searching for sequence ID matches!	
	==================================
	''')
	
	for filename in os.listdir(args.fasta_dir):
		if filename.endswith(".fasta") or filename.endswith(".fa"):
			fasta = list(SeqIO.parse(args.fasta_dir + "/" + filename, "fasta"))		
			for record in fasta:
				if record.description in seq_ids:
					ortho_seqs.append(record)
					
	# create the output fasta files				
	print("Matches found! Now writing orthogroup " + ortho_id + " nucleotide sequences as a fasta file")
	SeqIO.write(ortho_seqs, args.out_dir + "/" + ortho_id + ".fasta", "fasta")

print("Finished processing Orthofinder output!")
