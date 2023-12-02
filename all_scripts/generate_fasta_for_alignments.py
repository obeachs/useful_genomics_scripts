import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import sys
import os
from os import listdir
from Bio import pairwise2
import random


ATGGATAACACAAACCGTCTTCGCCGCCTTCACTGTCATAAACAACCCAAGTTCACTCATAGCTCTCAAGAAGTGAGTAGTATGAAATGGGAGTTTATCAATATGACCGAACAAGAAGAAGATCTCATCTTTAGAATGTACAGACTTGTTGGCGACAGGTGGGATTTAATAGCAAGAAGAGTGGTGGGACGTGAGGCAAAGGAGATAGAGAGATACTGGATTATGAGAAATTGTGACTATTTCTCCCACAAATGA
ATTCGTCCTTTCCTTCTGTATGTCACATATTTTGTTGTGTATGTCAACAACTAGGCTTTTTTGTTTGGTATCTATAATTCGTTCCTTTTTTTAGTTTGCCGGATCAGAATGTATTTTATGATTCATATCTTCTACAAAAGAGCCGGATATACTATATAGTTTCCTAGTTTATGTTAAAATTTATGTATCCACTTTGTTCCATCAAAGTTATTAAAAAGATAAAAGTTGCAAGATAAGAAAAAATAACTACATACAAATTGCAAACAAAAAAATAATATTTTATAATAAATTCCTCATTTTGAATACATAATAATATATCGTTTATAAATACAACATATATCATATGCATGTATTATTATTACAAAGTGGAATGATATCCATAATTTAAAATGTTTGTATTTGTGATCTGCTTTGCTTCTGAAAATTAGTTTTAGCATTTAATTTATTCAGTAAGTTATGGTGTACCTATATTTTGTAAGTTAAAGGTTATGTGGCCAAGTAGGAGAGTTGGGGACCTCTCTTTGTGTCTCTACATAGTTTTGTAAACCTCACATCTCTTTCTTCGCTTGCATTCTCCAAACTCCAATTTTTTTTGTTTCTCTCAATAATATTTGTGTTCATACTGTTTCGCTGTTTCCAATACATACTAATGGATAACACAAACCGTCTTCGCCGCCTTCACTGTCATAAACAACCCAAGTTCACTCATAGCTCTCAAGGTTTGTCTTTGATTTAGTTGGAAATCTAATTGTTCTTCCTTTTGGTCTACGCTCATTCATTTGACCTGTAACATTAAATTGTTTTGTGATGTACGTTTTCACTGATATTTTGAATTTATATACTAACTAAAATGCAGAAGTGAGTAGTATGAAATGGGAGTTTATCAATATGACCGAACAAGAAGAAGATCTCATCTTTAGAATGTACAGACTTGTTGGCGACAGGCACGTAACATTTTTCTTCTTCTGTTTATATATTGATCTAGTTCCAAATAGAAAAATTAAATTCTTGACCTAATTGTATTATAATATGTGTGAAAATATGTGCTTTAACAATTAAATATTTTTGTTTGTGTGTGAGGGATAAAACTGTTTGCCGTGATTCAAGCTTCTATCCAAGATCATTTTTACCTTGTTCTGAAAATTCGGATGTCGTGAGAGTCAGGATCTGTCCACACCAAATTTTGGATCATATATTATAATTGTTTTTATAAAAAAAAATATCCGGATCTATAGTTATATATATATTTTGAAAGTTTTTTTTACATAAAAATAATTTTAGATACCAAATTTATTTTAAAAACTACTTTTTAAACTTAATCCATATTTTCTAAGACTGTTCAATATTTAATTTTTATTTTTTAATAATTCGAATCTGAAGCTTATATAAAGAATCGATGACTTGCCTAACAATTGTTATAAACTATACTATAATTTTTTTTATTTTATAATACTTCCCTTTGATTCTAATTAATTCAAGTCTTACAAAAAATAAATTGTTTCATAAATATATCAAGTTTTCTATGCATTCTTTAGTTAATCATGAGATGAACTGTAAACTCGAAAAAAAATAATACTATTGGTTTAGAGTTAGAGAAATTTATTAATCACAAAATGATATATTTATAATCAAATTTTAATTATTCGTTCTTAATATGTGTTAAAACCTTAATACTTAAATTAATTAGAAATAGAAAAACTATTGTTGTTTTAAATTTATATAAATATATAGACTGTTAATGTAATGAAATCAAAAATGTAATTTTAGTTTTAATAATTCATTTTTTTACAATAATCAGAATACTAATAATAACATATATATGTTGGTCACAAATTGAATCATATGATCAGGTGGGATTTAATAGCAAGAAGAGTGGTGGGACGTGAGGCAAAGGAGATAGAGAGATACTGGATTATGAGAAATTGTGACTATTTCTCCCACAAATGACGCAAAAGTCACTCTCTTTAAAAAATCTCGAGAATCTCCCTGTTTTTCTATTTCTCATCCTCATTAAATAATTTCTGTGATTGTGGACAAGAAAAGCTTTATTAAGTTGTCGAGTGGGGTTCTCCTCATCATTAAATCTGTTATGGT

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--fasta', required = True, default = None, type = str)
parser.add_argument('--align', required = False, default=0,type= int)
args =parser.parse_args()


id_list = []
seq_list = []


fasta = args.fasta


with open(fasta,'r') as fa:
    for seq in SeqIO.parse(fa,'fasta'):
        id_list.append(str(seq.id))
        seq_list.append(str(seq.seq))
# function to get unique values
def unique(list1):
  
    # initialize a null list
    unique_list = []
      
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list


uni = unique(id_list)
newuni = []
for check in uni:
    if len(check) >= 9 and check[0:1] =='A' and 'cdna' not in check:
        newuni.append(check)
print(newuni)
'''This works well - finds the ID you're looking for then will scan forward and
take any other ids after that (the matches) until it reaches the next Arabidopsis
ID at which point it stops'''

for check in newuni:
    tair_id = check[0:9]
    newseq = []
    newid = []
    newseq_write = []
    newid_write = []
    print(check)
    location_list = [i for i,x in enumerate(id_list) if x[0:9]==tair_id]
    for i in location_list:
        newid.append(id_list[i])
        if 'cds' in id_list[i]:
            newid.append('SnAlba_' + tair_id + '_' +str(i))
        newseq.append(seq_list[i])
        if 'cds' in id_list[i]:
            newseq.append(seq_list[i+1])
    for i in range(len(newseq)):
        if newid[i] not in newid_write:
            newid_write.append(newid[i])
            newseq_write.append(newseq[i])
        else:
            if newseq[i] not in newseq_write:
                replace = newid[i] + '_' + str(i)
                newid_write.append(replace)
                newseq_write.append(newseq[i])

        


    namefasta = os.path.basename(fasta)
    result_fasta = namefasta.split('.')[0] +'_'+ tair_id + '_matches.fa'
    if len(newid_write) != len(newseq_write):
        print('ERROR')
        for id in newid_write:
            print(id)
        for seq in newseq_write:
            print(seq)
     
    with open(result_fasta,'w+') as result:
        for i in range(len(newid_write)):
            result.write('>' + newid_write[i] + '\n' + newseq_write[i] + '\n')
'''
if args.align ==1:
    aligned_fasta = fasta.split('.')[0] + '_' + check + '_aligned'
    align = AlignIO.read(result_fasta, "fasta")
    with open(aligned_fasta,'w+') as align_final:
        SeqIO.write(align, align_final, 'nexus')

'''
