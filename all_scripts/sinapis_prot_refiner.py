import Bio
import sys
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio.Seq import Seq



parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--fasta', required = True, default = None, type= str)
args = parser.parse_args()
def levenshtein_ratio_and_distance(s, t, ratio_calc = False):
    """ levenshtein_ratio_and_distance:
        Calculates levenshtein distance between two strings.
        If ratio_calc = True, the function computes the
        levenshtein distance ratio of similarity between two strings
        For all i and j, distance[i,j] will contain the Levenshtein
        distance between the first i characters of s and the
        first j characters of t
    """
    # Initialize matrix of zeros
    rows = len(s)+1
    cols = len(t)+1
    distance = np.zeros((rows,cols),dtype = int)

    # Populate matrix of zeros with the indeces of each character of both strings
    for i in range(1, rows):
        for k in range(1,cols):
            distance[i][0] = i
            distance[0][k] = k

    # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0 # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
            else:
                # In order to align the results with those of the Python Levenshtein package, if we choose to calculate the ratio
                # the cost of a substitution is 2. If we calculate just distance, then the cost of a substitution is 1.
                if ratio_calc == True:
                    cost = 2
                else:
                    cost = 1
            distance[row][col] = min(distance[row-1][col] + 1,      # Cost of deletions
                                 distance[row][col-1] + 1,          # Cost of insertions
                                 distance[row-1][col-1] + cost)     # Cost of substitutions
    if ratio_calc == True:
        # Computation of the Levenshtein Distance Ratio
        Ratio = ((len(s)+len(t)) - distance[row][col]) / (len(s)+len(t))
        return Ratio
    else:
         print(distance) # Uncomment if you want to see the matrix showing how the algorithm computes the cost of deletions,
        # insertions and/or substitutions

fasta = args.fasta

orig_out = fasta.split('.')[0] + '_refined.fa'
with open(fasta,'r') as fa, open(orig_out,'w+') as outfile:
    check_len = 0
    new_len = 0
    for seq in SeqIO.parse(fa,'fasta'):
        if seq.id[0] == 'A':
            if 'cdna' in seq.id:
                check_len = len(str(seq.seq))
                outfile.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
        else:
            new_len = len(str(seq.seq))
        if abs(new_len - check_len) > 12:
            continue
        else:
            outfile.write('>' + seq.id + '\n' + str(seq.seq) + '\n')

# with open(orig_out,'r') as orig:
#     tair_seqs = []
#     seen = []
#     for seq in SeqIO.parse(orig,'fasta'):
#         if seq.id[0] == 'A':
#             if 'cdna' in seq.id:
#                 tair_seqs.append(seq.id)

# print(tair_seqs)
# with open(orig_out,'r') as orig:
#     seen = []
#     for seq in SeqIO.parse(orig,'fasta'):
#         if tair_seqs.count(seq.id[0:8])>1:
#             print('CHECK A PASSED')
#             if 'cdna' in seq.id:
#                 if seq.id not in seen:
#                     print('CHECK B PASSED')
#                     new_main_id = str(seq.id)
#                     new_outname = orig_out.split('.')[0] + '_' + seq.id + '_refined.fa'
#                     seen.append[seq.id]
#                     with open(new_outname,'w+') as new_out:
#                         for seq in SeqIO.parse(orig,'fasta'):
#                             if seq.id == new_main_id:
#                                 print('CHECK C PASSED')
#                                 new_out.write('>' + new_main_id + '\n' + str(seq.seq) + '\n')
#                             if seq.id[0] != 'A':
#                                 if 'cdna' not in seq.id:
#                                     print('CHECK D PASSED')
#                                     new_out.write('>' + seq.id + '\n' + str(seq.seq) + '\n')

        

