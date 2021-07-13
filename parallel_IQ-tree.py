#!/usr/bin/env python

import multiprocessing
import glob
import os
import subprocess
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Parallelized gene tree inference using IQ-TREE")

parser.add_argument('--aln_dir', required=True, help='Directory containing aligned nucleotide sequence files in fasta format', action='store')
parser.add_argument('--out_dir', required=True, help='Output directory where IQ-TREE results will be written', action='store')
parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')

args=parser.parse_args()

#make the output directory if it doesn't already exist

if args.out_dir.endswith('/'):
    out_dir = args.out_dir
else:
    out_dir = args.out_dir + '/'
print(out_dir)

if not os.path.exists(out_dir):
        os.makedirs(out_dir)

alns = glob.glob(args.aln_dir + '*nt_ali.fasta')
print(alns)

def run_phylo(aln_file):
    # Run IQ-Tree ML phylo inference with model selection, ultrafast boostrap, and SH-aLRT test
    ID = aln_file.split('.')[1].strip('/')
    outprefix = out_dir + ID
    print(ID)
    print(outprefix) 
    #command1 = 'iqtree -s ' + aln_file + ' -m MFP -bb 1000 -alrt 1000 -pre ' + args.out_dir + '/' + ID
    #print(command1)
    subprocess.run(args=['iqtree', '-s', aln_file, '-m', 'MFP', '-nt', 'AUTO', '-ntmax', '5',
    '-bb', '1000', '-alrt', '1000', '-pre', outprefix])#, '-pre', out, ID'])


# Parallelize gene tree inference using multiprocessing

p = multiprocessing.Pool(processes=args.threads)

for aln in alns:
        print("Inferring gene tree for " + aln)
        #run_phylo(aln)
        p.apply_async(run_phylo, [aln])

p.close()
p.join() # Wait for all child processes to close.
