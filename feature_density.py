#!/usr/bin/env python

# use pybedtools to get gene (or other feature) density in sliding windows across a particular scaffold or chromosome

import pybedtools as pbt
from Bio import SeqIO
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Takes genome fasta and bed file of features and calculates number of features in sliding windows")

parser.add_argument('--fasta', required=True, help='genome assembly fasta file', action='store')
parser.add_argument('--genes', required=True, help='gff file of features genes', action='store')
parser.add_argument('--repeats', required=True, help='gff file of features repeats', action='store')
parser.add_argument('--window', required=True, help='size of windows to count features', action='store')
parser.add_argument('--step', required=True, help='size of the interval between windows', action='store')
parser.add_argument('--chrom', required=True, help='chromosome or scaffold to be used', action='store')

args=parser.parse_args()

# create variables needed to run the script from the arguments

genome = args.fasta
mygenes = args.genes
myrepeats = args.repeats
window_size = args.window
step_size = args.step
chrom = args.chrom

# specify output bed file names
gene_output = f"{chrom}_gene_density_{int(window_size) / 1000}kbwindow_{int(step_size) / 1000}kbstep.bed"
repeat_output = f"{chrom}_repeat_density_{int(window_size) / 1000}kbwindow_{int(step_size) / 1000}kbstep.bed"

# create bed file from fasta genome

genome_list = []

for record in SeqIO.parse(genome, 'fasta'):
    header = record.id
    length = str(len(record.seq) - 1)
    start = str(0)
    line = '\t'.join((header, start, length))
    genome_list.append(line)

# make bed file and filter results to include only the desired scaffold
genome_bed = pbt.BedTool(genome_list).filter(lambda x: x[0] == chrom).saveas()

genes_bed = pbt.BedTool(mygenes).filter(lambda x: x[0] == chrom).saveas()
repeats_bed = pbt.BedTool(myrepeats).filter(lambda x: x[0] == chrom).saveas()

# get sliding windows for the desired scaffold given window and step size

windows = pbt.BedTool().window_maker(s=step_size, w=window_size, b=genome_bed).saveas()

# use bedtools coverage function to get the proportion (7th column) of each window covered by the feature
# require that a window has to be at least half of the window size to remove artifacts at the end of the scaffold
gene_coverage = windows.coverage(b=genes_bed).filter(lambda x: int(x[5]) >= int(window_size) / 2).saveas(gene_output)

repeat_coverage = windows.coverage(b=repeats_bed).filter(lambda x: int(x[5]) >= int(window_size) / 2).saveas(repeat_output)

# read results bed files to make a dataframe for plotting
gene_df = pd.read_csv(gene_output, sep = '\t', header = None, names = ['Chrom', 'Start', 'Stop', 'Overlapping features', 'Covered bases', 'Window size', 'Gene fraction coverage'])

repeat_df = pd.read_csv(repeat_output, sep = '\t', header = None, names = ['Chrom', 'Start', 'Stop', 'Overlapping features', 'Covered bases', 'Window size', 'Repeat fraction coverage'])

def line_plot(gene_df, repeat_df):
	sns.set(rc={'figure.figsize':(8,2)})
	sns.set_style("white")
	sns.set_context("paper", font_scale=0.75)

	ax1 = sns.lineplot(data = gene_df, x = 'Start', y = 'Gene fraction coverage', color = 'cornflowerblue')
	ax2 = sns.lineplot(data = repeat_df, x = 'Start', y = 'Repeat fraction coverage', color = 'darkorchid')
	plt.xlabel('Scaffold Position')
	plt.ylabel("Genomic element density")
	plt.legend(['genes', 'repeats'], loc = 'upper right')
	plt.savefig(f"{chrom}_feature_density.png", bbox_inches='tight')
	plt.savefig(f"{chrom}_feature_density.eps", bbox_inches='tight')
	plt.show()
	plt.clf()

line_plot(gene_df, repeat_df)
