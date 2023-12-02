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



# tab = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/joi_assembly_info.txt')
# joi_count = tab['Reads'].to_list()
# read_len = tab['Cutoff'].to_list()
# with open('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/joi_readlength_counts.txt', 'w+') as oop:
#     oop.write('Length' + '\t' + 'Count' + '\n')
#     for i in range(len(joi_count)-1):
#         oop.write(str(read_len[i]) + '\t' + str((joi_count[i]) - (joi_count[i+1])) + '\n')
# '''




def fasta_condenser(fasta, tair=0):
    #Preferably use this on the fasta file that has less info/annotation
    #with it - usually the query fasta
    namelist = []
    seqlist = []
    if tair==0:
        with open (fasta,'r') as fa:
            for seq in SeqIO.parse(fa,'fasta'):
                namelist.append(seq.id)
                seqlist.append(str(seq.seq))
        df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
    else:
        with open (fasta,'r') as fa:
            for seq in SeqIO.parse(fa,'fasta'):
                namelist.append(seq.id[0:9])
                seqlist.append(str(seq.seq))
        df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
    return df

fasta = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/all_nanopores_combind5000_k17k25k31.fa'
reads = fasta_condenser(fasta)
print(str(reads['Seq'].min()))
seqs = reads['Seq'].to_list()
l = []
for i in seqs:
    l.append(len(i))
print(l)
maximum = max(l)
maximum = round(maximum,-2)
print(maximum)
num = int(maximum/1000)
print(num)

output = []
output_two = []
for i in range(num):
    x = i * 1000
    y = (i-1) * 1000
    print(x)
    print(y)
    output.append(sum(y <= float(n) <= x for n in l))
    output_two.append(x)
with open('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/readlength_counts_1000.txt', 'w+') as oop:
    oop.write('Length' + '\t' + 'Count' + '\n')
    for i, j in zip(output, output_two):
        oop.write(str(j) + '\t' + str(i) + '\n')



