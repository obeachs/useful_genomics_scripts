import numpy as np
import pandas as pd
import os
import sys
import argparse
import Bio
from Bio import Align
from Bio.Seq import Seq

from Bio import SeqIO

#Seeting up the parameters for the alignment, nothing fancy - just want to make sure it's
#local so that the score isn't always taking the full length of the sequence into account
#Using +1 and -1 for each of the scores seems to be okay.

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--f', required = True, type = str)
parser.add_argument('--r', required = True, type = str)
parser.add_argument('--fasta', required =True, type=str)
args =parser.parse_args()


def fasta_condenser(fasta, tair=0):
    '''Preferably use this on the fasta file that has less info/annotation
    with it - usually the query fasta'''
    namelist = []
    seqlist = []
    if tair==0:
        with open (fasta,'r') as fa:
            for seq in SeqIO.parse(fa,'fasta'):
                namelist.append(seq.description)
                seqlist.append(str(seq.seq))
        df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
    else:
        with open (fasta,'r') as fa:
            for seq in SeqIO.parse(fa,'fasta'):
                namelist.append(seq.id[0:9])
                seqlist.append(str(seq.seq))
        df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
    return df


fasta = fasta_condenser(args.fasta)
f = args.f
r = args.r


#For each dna sequence in the 'Seq' column of fasta, check if the i string is present in it, allowi for up to 4 mismatches
#If it is present, add the TAIR_ID to the list 'matches'
matches = []
#If the score is greater than the length of the sequence -4, then it is likely a match
faim = len(f)
raim = len(r)

for j in range(len(fasta['Seq'])):
    falignments= aligner.align(fasta['Seq'][j], f)
    ralignments= aligner.align(Seq(fasta['Seq'][j]).reverse_complement(),r)
    if falignments.score > faim-5 and ralignments.score > raim-5:
        print(fasta['ID'][j])
        print('Forward primer:  ' + str(falignments.score))
        print('Reverse primer:  ' + str(ralignments.score))
        print(falignments[0])
        print(ralignments[0])
