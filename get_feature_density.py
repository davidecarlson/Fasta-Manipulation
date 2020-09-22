#!/usr/bin/env python

# use pybedtools to get gene (or other feature) density in sliding windows across a particular scaffold or chromosome

import pybedtools as pbt
from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser(description="Takes genome fasta and bed file of features and calculates number of features in sliding windows")

parser.add_argument('--fasta', required=True, help='genome assembly fasta file', action='store')
parser.add_argument('--bed', required=True, help='bed file of features (e.g., genes)', action='store')
parser.add_argument('--window', required=True, help='size of windows to count features', action='store')
parser.add_argument('--step', required=True, help='size of the interval between windows', action='store')
parser.add_argument('--chrom', required=True, help='chromosome or scaffold to be used', action='store')
parser.add_argument('--output', required=True, help='name of output bed file to write', action='store')

args=parser.parse_args()

# create variables needed to run the script from the arguments

genome = args.fasta
mybed = args.bed
window_size = args.window
step_size = args.step
outfile = args.output
chrom = args.chrom

# create bed file from fasta genome

genome_list = []

for record in SeqIO.parse(genome, 'fasta'):
    header = record.id
    length = str(len(record.seq) - 1)
    start = str(0)
    line = '\t'.join((header, start, length))
    genome_list.append(line)

# make bed file and filter results to include only the desired scaffold
genome_bed = pbt.BedTool(genome_list).filter(lambda x: x[0] == chrom)
#print(genome_bed)

# subset the features bed file to only include the desired scaffold
# more efficient to do this with the .filter method in pybedtools instead of with regular expressions

#with open(mybed) as file:
#	myscaffold = [line for line in file if re.match('^HiC_scaffold_16\t', line) is not None]
#scaffold_bed = pbt.BedTool(myscaffold)

feature_bed = pbt.BedTool(mybed).filter(lambda x: x[0] == chrom)

# get sliding windows for the desired scaffold

windows = pbt.BedTool().window_maker(s=step_size, w=window_size, b=genome_bed)

# count number of genes in each window with intersect method and write results to new bed file
# clip off chrom end where the region drops below the window size

#gene_count = windows.intersect(b=feature_bed, c=True).filter(lambda x: int(x[2]) - int(x[1]) == int(window_size)).saveas(outfile)
#print(gene_count)

gene_count = windows.intersect(b=feature_bed, c=True).saveas(outfile)
