#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Take fasta file as input and print the lengths of each sequence')

parser.add_argument('--fasta', required=True, dest='fasta', action='store')

args=parser.parse_args()

for record in SeqIO.parse(args.fasta, "fasta"):
	print(record.id + "\t" + str(len(record.seq)))
