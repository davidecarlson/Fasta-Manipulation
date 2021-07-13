import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Take fastq file as input and convert to fasta')

parser.add_argument('--fastq', required=True, dest='fastq', action='store')
parser.add_argument('--fasta', required=True, dest='fasta', action='store')

args=parser.parse_args()

fastq = args.fastq
fasta = args.fasta

count = SeqIO.convert(fastq, "fastq", fasta, "fasta")
print(f"Converted {count} records")
