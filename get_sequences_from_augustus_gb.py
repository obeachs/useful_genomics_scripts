import re
import sys
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
import Bio
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir



# outfile = '/Volumes/sesame/joerecovery/genomes/tair_id_list.txt'
# gff = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/TAIR10_GFF3_genes.txt')
# gff2 = gff.loc[gff['TYPE'] == 'gene']
# list = gff['ID2'].to_list()

# newlist = []
# for i in range(len(list)):
#     start = list[i]
#     startpoint = start.find('=')
#     newlist.append(list[i][startpoint + 1: startpoint+10])


# newerlist = []
# for l in newlist:
#     if l not in newerlist:
#         newerlist.append(l)

# with open(outfile,'w+') as out:
#     for i in newerlist:
#         out.write(i + '\n')







# def fasta_condenser(fasta, tair=0):
#     '''Preferably use this on the fasta file that has less info/annotation
#     with it - usually the query fasta'''
#     namelist = []
#     seqlist = []
#     if tair==0:
#         with open (fasta,'r') as fa:
#             for seq in SeqIO.parse(fa,'fasta'):
#                 namelist.append(seq.id)
#                 seqlist.append(str(seq.seq))
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
#     else:
#         with open (fasta,'r') as fa:
#             for seq in SeqIO.parse(fa,'fasta'):
#                 namelist.append(seq.id[0:9])
#                 seqlist.append(str(seq.seq))
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
#     return df

# tb = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/swissprot_blastx/all_rna_seq_merged_maingenome_dedup.txt')
# ref = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/all_rna_seq_merged/all_rna_seq_merged_maingenome.fa')
# all_out = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/swissprot_blastx/all_rna_seq_merged_maingenome_all_swissprot_gene_hits.fa'
# tair_out_no_dups = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/swissprot_blastx/all_rna_seq_merged_maingenome_tair_swissprot_gene_hits_no_dups.fa'
# print(ref)

# df = ref.merge(tb, left_on='ID', right_on='sinapis_id', how='inner',copy=False)
# print(df)
# df = df[['sinapis_id','tair_id','Seq']]
#                 print(seq.seq[start:end])
# print(df)
# tair_df = df[df['tair_id'].str.contains('ARATH')]
# tair_df = tair_df.drop_duplicates(subset='sinapis_id')
# print(tair_df)
# tair = tair_df['tair_id'].to_list()
# sin = tair_df['sinapis_id'].to_list()
# seq = tair_df['Seq'].to_list()
# with open(tair_out_no_dups,'w+') as outy:
#     for i in range(len(seq)):
#         id =sin[i]
#         outy.write('>' + id + '\n' + seq[i] + '\n')


# with open(ref,'r') as reference, open(all_out,'w+') as outy:
#     for i in tb['sinapis_id']:
#         for seq in SeqIO.parse(reference,'fasta'):
#             if seq.id == i:
#                 print(i)


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--id', required = True, type = str)
parser.add_argument('--start',required = False, type = int, default=0)
parser.add_argument('--end', required = False, type= int, default=0)
parser.add_argument('--ref', required = False)
parser.add_argument('--rev', required = False, type= int, default=0)
args = parser.parse_args()

id = args.id
start = args.start
ref = args.ref
end = args.end
if end == 0:
    end = -1 
rev = args.rev
# end = end+500
with open(ref,'r') as reference:
    for seq in SeqIO.parse(reference,'fasta'):
        if seq.id == id:
            print(seq.id)
            print(seq.seq[start:end])
