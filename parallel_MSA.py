#!/usr/bin/env python

import multiprocessing
import glob
import os
import subprocess
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Parallelized alignment of multiple fasta files using TranslatorX")

parser.add_argument('--fasta_dir', required=True, help='Directory containing nucleotide coding sequence files in fasta format', action='store')
parser.add_argument('--out_dir', required=True, help='Output directory where orthogroup fasta files will be written', action='store')
parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')

args=parser.parse_args()

#make the output directory if it doesn't already exist


if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

fastas = glob.glob(args.fasta_dir + '*.fasta')
print(fastas)

def run_msa(fa_file):
    ID = fa_file.split('.')[1].strip('/')
    command1 = 'translatorx_vLocal.pl -i ' + fa_file + ' -p M -o ' + args.out_dir + '/' + ID
    #print(command1)
    subprocess.run(command1, stdout = True, shell = True)

# Parallelize the alignment of each fasta file multiprocessing

p = multiprocessing.Pool(processes=args.threads)

for fasta in fastas:
    num_seqs = len(list(SeqIO.parse(fasta, "fasta")))
    # ignore fasta file if it has less than 10 sequences
    if num_seqs >= 10:
        
        #run_msa(fasta)
        print("aligning " + fasta + " with " + str(num_seqs) + " sequences!")
        p.apply_async(run_msa, [fasta])
    else:
        print("Too few sequences in Orthogroup" + fasta)
    

p.close()
p.join() # Wait for all child processes to close.
