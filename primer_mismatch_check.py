import numpy as np
import pandas as pd
import os
import sys
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import Bio
from Bio import Align
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
parser.add_argument('--i', required = True, type = str)
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
i = args.i


#For each dna sequence in the 'Seq' column of fasta, check if the i string is present in it, allowi for up to 4 mismatches
#If it is present, add the TAIR_ID to the list 'matches'
matches = []
#If the score is greater than the length of the sequence -4, then it is likely a match
aim = len(i)
print('Primer length is ' + str(aim) + '\n')
for j in range(len(fasta['Seq'])):
    alignments = alignments = aligner.align(fasta['Seq'][j], i)
    if alignments.score > aim-4:
        print(fasta['ID'][j])
        print(alignments[0])
        print(alignments.score)
