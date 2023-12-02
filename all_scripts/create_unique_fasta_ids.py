import argparse
import Bio
import numpy
from Bio import SeqIO




parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--fasta', required = True)
args =parser.parse_args()

fasta = args.fasta

outfile = fasta.split('.')[0] + '_edited.fa'

with open(fasta,'r') as file, open(outfile,'w+') as out:
    seen = []
    count =0 
    for seq in SeqIO.parse(file,'fasta'):
        count += 1
        if str(seq.id) not in seen:
            seen.append(str(seq.id))
            out.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')
        else:
            out.write('>' + str(seq.id) + '_' + str(count) + '\n' + str(seq.seq) + '\n')

