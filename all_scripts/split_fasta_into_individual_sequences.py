import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
import argparse
import pysam

'''Splits a given fasta file containing N number of reads into N number of files
- Very useful when doing BLAST anaylses or need to upload individual sequences 
to sites
Usage split_fasta_individuals.py --fasta FASTAFILE --directory Folde where you want the sequences stored
'''
parser = argparse.ArgumentParser(description = 'Split fasta file into reads\n')
parser.add_argument('--fa', required = True, type=str)
parser.add_argument('--dir', required = True, type=str)

args =parser.parse_args()
fa = args.fa
folder = args.dir

fastafile = open(fasta,'r')
for seq in SeqIO.parse(fastafile,'fasta'):
    file = open(folder + '/' + str(seq.id) + '.fa', 'w+'):
    file.write('>' + seq.id + '\n' + str(seq.seq))
