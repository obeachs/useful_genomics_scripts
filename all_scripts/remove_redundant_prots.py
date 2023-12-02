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
import random
from Bio import pairwise2 as pw2

def df_to_write_fasta(df, outfile):
    with open(outfile,'w+') as outy:
        names = list[df.columns]
        l = len(df['ID'])
        print(l)
        for i in range(l):
            outy.write('>' + df.iloc[i]['ID'] + '\n' + df.iloc[i]['Seq'] +'\n')

def sequence_similarity(first_seq, second_seq):
    global_align = pw2.align.globalxx(first_seq, second_seq)

    seq_length = min(len(first_seq), len(second_seq))
    matches = global_align[0][2]

    percent_match = (matches / seq_length) * 100
    return(percent_match)


def fasta_condenser(fasta, tair=0):
    '''Preferably use this on the fasta file that has less info/annotation
    with it - usually the query fasta'''
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
    

protein_fasta = '/Volumes/sesame/joerecovery/genomes/brapa/protein.faa'
prot = fasta_condenser(protein_fasta)
#Very easy way to extract random reads from the condensed fasta dataframe
num = 10
new_prot = prot.sample(num)



sequences = new_prot['Seq'].to_list()
ids = new_prot['ID'].to_list()
out = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/gene_prediction/augustus_training/non_redundant_protein.faa'

remove_me = []
for seq in sequences:
    for comp_seq in sequences:
        if seq == comp_seq:
            continue
        if sequence_similarity(seq, comp_seq) >= 70:
            remove_me.append(seq)
            print('Found a sequence to remove: ' + seq)
            break
    continue
print(remove_me)

non_redundant_prot = new_prot[~new_prot['Seq'].isin(remove_me)]
print(non_redundant_prot)
print(non_redundant_prot.iloc[2]['Seq'])
df_to_write_fasta(non_redundant_prot, out)

#new_ids = ids
#new_seq = sequences

# for a, b in itertools.combinations(sequences, 2):
#     if sequence_similarity(a,b) > 70:
#         num = random.randint(0, 1)
#         if num == 1:
#             print(a)
#             location = sequences.index(a)
#             #new_seq.remove(a)
#             seq_del.append(a)
#             id_del.append(ids[location])
#         else:
#             print(b[1:-2])
#             location = sequences.index(b)
#             seq_del.append(b)
#             #new_seq.remove(b)
#             id_del.append(ids[location])


# for i in sequences[:]:
#    if i in seq_del:
#       sequences.remove(i)
# for j in ids[:]:
#    if j in id_del:
#       ids.remove(j)


# with open(out,'w+') as outfile:
#    for i in range(len(new_ids)):
#       outfile.write('>' + ids[i] + '\n' + sequences[i] + '\n')