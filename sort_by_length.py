#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO

# lightly adapted from the biopython tutorial found at http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc287

parser = argparse.ArgumentParser(description='Sort input fasta file by sequence length, from longest to shortest')

parser.add_argument('--fasta', required=True, help='input fasta file', action='store')

args=parser.parse_args()

# Don't want to read entire file into memory, so first just get the length and id of each sequence as a list and sort by length

len_and_ids = sorted(
    (len(rec), rec.id) for rec in SeqIO.parse(args.fasta, "fasta")
)

# reverse the order

ids = list(reversed([id for (length, id) in len_and_ids]))

#get each sequence record using the index function so we don't read the whole file into memory at once
record_index = SeqIO.index(args.fasta, "fasta")

records = (record_index[id] for id in ids)

# write sorted sequences to stdout

SeqIO.write(records, sys.stdout, "fasta")