import multiprocessing
import glob
import os
import argparse 
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Take Orthofinder output file and produce separate nucleotide fasta files for each orthogroup")

parser.add_argument('--orthogroups', required=True, help='Orthogroups.tsv file produced by Orthofinder', action='store')
parser.add_argument('--fasta_dir', required=True, help='Directory containing nucleotide coding sequence files in fasta format', action='store')
parser.add_argument('--out_dir', required=True, help='Output directory where orthogroup fasta files will be written', action='store')
parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')

args=parser.parse_args()

def process_orthofinder(row):

    ortho_id = row.split('\t')[0]
    seq_ids = row.split(ortho_id)[1].replace('\t',',')[1:].replace(' ', '').rstrip().split(',')

    # create empty list of sequences
    ortho_seqs = []

    # For each record from the input fasta files, check if record description matches an element in the 
    # sequence IDs list, and if so keep it #

    print(
    '''
    =================================
    Searching for sequence ID matches!
    ==================================
    ''')

    for record in fasta_records:
        if record.description in seq_ids:
            ortho_seqs.append(record)
            
    # if matches are not found, issue a warning
    if len(ortho_seqs) == 0:
        print("Warning! No matches found for " + ortho_id + "!")

    else:
        # create the output fasta files
        print("Matches found! Now writing orthogroup " + ortho_id + " nucleotide sequences as a fasta file")
        SeqIO.write(ortho_seqs, args.out_dir + "/" + ortho_id + ".fasta", "fasta")

#make the output directory if it doesn't already exist


if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

# read all fasta files into memory to prevent having to loop through them over and over

print("Reading input fastas.  This may take a moment!")

# create empty list of fasta records, which we will populate after reading each fasta file

fasta_records = []

for filename in os.listdir(args.fasta_dir):
    if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".cds"):
        fasta = list(SeqIO.parse(args.fasta_dir + "/" + filename, "fasta"))
        for record in fasta:
            fasta_records.append(record)

# quit if we don't have any recognizable input fasta files

if len(fasta_records) == 0:
    sys.exit("Error! No suitable fasta files found!")

# open the Orthogroups file

file =  open(args.orthogroups, 'r')
next(file, None) # skip header

# Parallelize the search of each orthogroup in the input file using multiprocessing
    
p = multiprocessing.Pool(processes=args.threads)
for row in file:    # launch a process for each row (up to the number of threads specified)
    p.apply_async(process_orthofinder, [row]) 
    

p.close()
p.join() # Wait for all child processes to close.
