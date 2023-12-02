import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2
import random





string = 'BINGUS'
print(string.split('/')[-1])
# ids = ['AT4G18960','AT5G53200','AT2G30432','AT2G46410','AT3G27920','AT1G79840','AT5G41315']
# df = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/all_nanopores_combind5000_k17k25k31_blast_to_TAIR10_cds_20101214_updated_dedup.txt')

# df['ids'] = df['tair_id'].str[0:9]
# print(df)
# df1 = df[df['ids'].isin(ids)]
# df1 = df1.drop(columns=['hit_length', 'tair_start', 'ids'])
# df1['Length'] = df1['sinapis_end'] - df1['sinapis_start']
# df1 = df1.sort_values('tair_id')
# df1.to_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/all_nanopores_combind5000_k17k25k31_blast_to_TAIR10_cds_20101214_updated_dedup_trichomes_only.txt', index=False)
# outfile ="/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/phylogeny/ATCG00120.1_same_length.fa"
# with open("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/phylogeny/ATCG00120.1.fa",'r') as fa, open(outfile,'w+') as out:
#     for seq in SeqIO.parse(fa,'fasta'):
#         list.append(len(seq.seq))
#         out.write('>' + str(seq.id) + '\n' + str(seq.seq[0:606]) + '\n')      

# print(min(list))


'''if 1 >2 andl
seque = Seq('tggatagacgacgtcggagacagagcaaggccaaagcgtcgtgttccgaagaagtgagtagcatcgaatgggaagctgtgaagatgactgaggaagaagaagatctcatttctcggatgtataaacttgtgggagacaggtgggaattgatagccggaaggataccgggacggacgccagaggagattgaaagatattggcttatgaaacatggtctcgtttttgccaacagacgaagagattttgttaggagatga')
print(seque.reverse_complement().complement())
nano = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/all_nanopores_combind5000_k17k25k31.fa'
outfile = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/all_nanopores_combind5000_k17k25k31_readlengths.txt'
with open(nano, 'r') as fasta, open(outfile,'w+') as out:
    for seq in SeqIO.parse(fasta,'fasta'):
        out.write(str(len(seq.seq)) + '\n')
       '''

'''
# Using readlines()
file1 = open('/Users/josephbeegan/Desktop/full_transposons.gff', 'r')
Lines = file1.readlines()


with open('/Users/josephbeegan/Desktop/clean_transposons.gff','w+') as out:
    out.write('chr' + '\t' + 'start' + '\t' + 'feature' + '\n')
    for line in Lines:
        start = line.split('\t')[0]
        chr = line.split('\t')[1]
        desc = line.split('\t')[2]
        id_start = desc.find('ID=')
        id  = desc[id_start + 3:id_start +12]
        out.write(str(start) + '\t' +chr + '\t' + id + '\n')






'''



# '''PLAN:

# Get blast matches(ids) -> extract matching ids -> put all matching ids into
# individual fasta files -> produce alignments'''

# table = pd.read_table('/Users/josephbeegan/Documents/sinapis_all_rnaseq_reads/All_blast_to_TAIR10_cds_20101214_updated_matches_list_dedup.txt')
# query_fasta = '/Users/josephbeegan/Documents/sinapis_all_rnaseq_reads/all_reads_transcript.fa'
# outfasta = '/Users/josephbeegan/Documents/sinapis_all_rnaseq_reads/matched_genes.fa'
# db_fasta = '/Users/josephbeegan/Documents/sinapis_temp/TAIR10_cds_20101214_updated.fasta'
# trichomes = '/Users/josephbeegan/Documents/sinapis_temp/trichome_ids.txt'


# def find_string(txt, str1):
#     return txt.find(str1, txt.find(str1)+1)


# '''Bit of a hacky way to get the symbol of the gene - lots of genes don't have
# symbols so hopefully leaving them as empty spaces is ok'''
# def get_symbol_from_tair(seq_description):
#     symbols_start = str(seq.description).find('Symbols:')+8
#     symbols_end = find_string(str(seq.description),'|')
#     symbol = str(seq.description)[symbols_start:symbols_end]
#     return symbol

# def fasta_condenser(fasta, tair=0):
#     '''Preferably use this on the fasta file that has less info/annotation
#     with it - usually the query fasta'''
#     namelist = []
#     seqlist = []
#     with open (fasta,'r') as fa:
#         for seq in SeqIO.parse(fa,'fasta'):
#             namelist.append(seq.id)
#             seqlist.append(str(seq.seq))
#     if tair==0:
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
#     else:
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
#     return df


# def extract_genelist(listfile):
#     with open(listfile, "r") as ids:
#         pre_ids = ids.read()
#         id_list = pre_ids.split("\n")


# query_df = fasta_condenser(query_fasta)
# db_df = fasta_condenser(db_fasta,tair=1)
# merged = db_df.merge(table, how='left',left_on='TAIR_ID', right_on='tair_id')
# merged = merged[["TAIR_ID","TAIR_Seq","sinapis_id"]]
# merged = merged.dropna()
# merged  = merged.merge(query_df,how='left',left_on = 'sinapis_id', right_on='ID')

# print(merged['TAIR_Seq'])

# final_tair_seq = merged['TAIR_Seq'].tolist()
# final_tair_id = merged['TAIR_ID'].tolist()
# final_query_id = merged['ID'].tolist()
# final_query_seq = merged['Seq'].tolist()


# with open(outfasta,'w+') as outfile:
# 	for i in range(len(final_tair_id)):
# 		outfile.write('>' + final_tair_id[i] + '\n' + final_tair_seq[i] + '\n' + '>' +
# 		final_query_id[i] + '\n' + final_query_seq[i] + '\n')






# # with open(db_fasta,'r') as db:
# #     for seq in SeqIO.parse(db,'fasta'):
# #         try:
# #             found_id = table.loc[table['tair_id'] == str(seq.id), 'sinapis_id'].values[0]
# #             found_seq = fasta_df.loc[table['ID'] == str(found_id), 'Seq'].values[0]
# #             print(str(seq.id_ +'\n' +str(seq.seq) +'\n' + found_id + '\n' + found_seq))
# #         except:
# #             continue



# #
# #
# #
# #
# # parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
# # parser.add_argument('--query', required = True, type = str)
# # parser.add_argument('--out', required = True, type = str)
# # args =parser.parse_args()
# # query = args.query
# # out = args.out
# # with open(query,'r') as xmlblast, open(out,'w+') as outfile:
# #     soaptair_record = NCBIXML.parse(xmlblast)
# #     for record in soaptair_record:
# #         for description in record.descriptions:
# #             for alignment in record.alignments:
# #                 for hsp in alignment.hsps:
# #                     tair_id =alignment.title.split(' ')[0]
# #                     barley_id = record.query
# #                     end = hsp.query_start + alignment.length
# #                     outfile.write(tair_id + '\t' + record.query + '\n')
# #
# #
# #
# # table = pd.read_table('/Users/josephbeegan/Documents/sinapis_temp/default_blast_analysis/sinapis_to_tair_protein_analysis_results_dedup_with_sequences.txt')
# # print(table['sinapis_id'])
# # id_list = table['sinapis_id'].to_list()
# # print(id_list)
# #
# # match_TCL1_seq = 'CAAGTCTGTACGTTCGAGAGATGAGATCTTCTTCTTGTTCGGTCATATTGATAAACTCCCATTTGATAATAGTCACTTCTGCATTTTAGTTAGTTATACGCAACCAATAATCAAATATATAAAATAATGAGGGAGAGAAGCTAATTTAGTAAAAAGAATTTAGGTAAAAATAAGTTCAGGGTAGAGGATGAGTATAGACCCAAAGGAAGAACCATTAAATATTAACGAAATCAAAGAGAGACCTTCGGAGTTATA'
# # scaffolds = open('/Users/josephbeegan/Documents/sinapis_temp/default_blast_analysis/sinapis_scaffolds_with_tair_matches.fa','r')
# # for seq in SeqIO.parse(scaffolds,'fasta'):
# #     if seq.id == 'scaffold_22':
# #         scaffold_22 = seq.seq
# # start_pos = int(scaffold_22.find(match_TCL1_seq)) - 100
# # end_pos = int(scaffold_22.find(match_TCL1_seq)) + len(match_TCL1_seq) + 100
# # potential = scaffold_22[start_pos:end_pos]
# # print(potential)
# # print(potential.find('ATG'))
# #
# #
# # #
# # # ks_seq=Seq('ATGGAGAGATCAAACAGCATAGAGTTGAGGAACAGCTTCTATGGCCGTGCAAGAACTTCACCATGGAGCTATGGAGATTATGATAATTGCCAACAGGATCATGATTATCTTCTAGGGTTTTCATGGCCACCAAGATCCTACACTTGCAGCTTCTGCAAAAGGGAATTCAGATCGGCTCAAGCACTTGGTGGCCACATGAATGTTCACAGAAGAGACAGAGCAAGACTCAGATTACAACAGTCTCCATCATCATCTTCAACACCTTCTCCTCCTTACCCTAACCCTAATTACTCTTACTCAACCATGGCAAACTCTCCTCCTCCTCATCATTCTCCTCTAACCCTATTTCCAACCCTTTCCTCCTCCATCCTCACCAAGATATAGGGCAGGTTTGATCCGTTCCTTGAGCCCCAAGTCAAAACATACACCAGAAAACGCTTGTAAGACTAAGAAATCATCTCTTTTAGTGGAGGCTGGAGAGGCTACAAGGTTCACCAGTAAAGATGCTTGCAAGATCCTGAGGAATGATGAAATCATCAGCTTGGAGCTTGAGATTGGTTTGATTAACGAATCAGAGCAAGATCTGGATCTAGAACTCCGTTTGGGTTTCGCTTAA')
# # #
# # # wt_seq=Seq('ATGGAGAGATCAAACAGCATAGAGTTGAGGAACAGCTTCTATGGCCGTGCAAGAACTTCACCATGGAGCTATGGAGATTATGATAATTGCCAACAGGATCATGATTATCTTCTAGGGTTTTCATGGCCACCAAGATCCTACACTTGCAGCTTCTGCAAAAGGGAATTCAGATCGGCTCAAGCACTTGGTGGCCACATGAATGTTCACAGAAGAGACAGAGCAAGACTCAGATTACAACAGTCTCCATCATCATCTTCAACACCTTCTCCTCCTTACCCTAACCCTAATTACTCTTACTCAACCATGGCAAACTCTCCTCCTCCTCATCATTCTCCTCTAACCCTATTTCCAACCCTTTCTCCTCCATCCTCACCAAGATATAGGGCAGGTTTGATCCGTTCCTTGAGCCCCAAGTCAAAACATACACCAGAAAACGCTTGTAAGACTAAGAAATCATCTCTTTTAGTGGAGGCTGGAGAGGCTACAAGGTTCACCAGTAAAGATGCTTGCAAGATCCTGAGGAATGATGAAATCATCAGCTTGGAGCTTGAGATTGGTTTGATTAACGAATCAGAGCAAGATCTGGATCTAGAACTCCGTTTGGGTTTCGCTTAA')
# # #
# # # sup_1 = Seq('ATGGAGAGATCAAACAGCATAGAGTTGAGGAACAGCTTCTATGGCCGTGCAAGAACTTCACCATGA')
# # #
# # # print(ks_seq.translate())
# # # print(sup_1.translate())
# # # print(wt_seq.translate())
# # # print(len('MERSNSIELRNSFYGRARTSPWSYGDYDNCQQDHDYLLGFSWPPRSYTCSFCKREFRSAQALGGHMNVHRRDRARLRLQQSPSSSSTPSPPYPNPNYSYSTMANSPPPHHSPLTLFPTLSSSILTKI'))
