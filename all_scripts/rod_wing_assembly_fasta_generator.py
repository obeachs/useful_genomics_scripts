import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2



tair_cds_fasta = '/Users/josephbeegan/Documents/sinapis_temp/TAIR10_cds_20101214_updated.fasta'
sinapis_fasta = '/Users/josephbeegan/Documents/sinapis_temp/sinapis_alba_var_s2_gc0560-79.maingenome.fasta'



with open(sinapis_fasta,'r') as sin:
    for seq in SeqIO.parse(sin,'fasta'):
        if seq.description == 'scaffold_22':
            print(seq.seq)





# sinapis_homologue_table = pd.read_table('/Users/josephbeegan/Documents/sinapis_temp/default_blast_analysis/All_blast_to_TAIR10_cds_20101214_updated_dedup.txt')
# print(sinapis_homologue_table)
# sinapis_homologue_table['tair_id'] = sinapis_homologue_table['tair_id'].str[0:9]
# '''Genes can be changed, in this case it is:
# AG, TRY, TCL1, CPC, GL1, GL2, GL3, SUP'''
# id_list =['AT4G18960','AT5G53200','AT2G30432','AT2G46410','AT3G27920','AT1G79840','AT5G41315','AT3G23130']
# wanted_table = sinapis_homologue_table[sinapis_homologue_table['tair_id'].isin(id_list)]
#
#
#
# def find_nth(string, character, n):
#     start = string.find(character)
#     while start >= 0 and n > 1:
#         start = string.find(character, start+len(character))
#         n -= 1
#     return start
#
# def getIndexes(dfObj, value):
#
#     # Empty list
#     listOfPos = []
#
#     # isin() method will return a dataframe with
#     # boolean values, True at the positions
#     # where element exists
#     result = dfObj.isin([value])
#
#     # any() method will return
#     # a boolean series
#     seriesObj = result.any()
#
#     # Get list of column names where
#     # element exists
#     columnNames = list(seriesObj[seriesObj == True].index)
#
#     # Iterate over the list of columns and
#     # extract the row index where element exists
#     for col in columnNames:
#         rows = list(result[col][result[col] == True].index)
#
#         for row in rows:
#             listOfPos.append((row, col))
#
#     # This list contains a list tuples with
#     # the index of element in the dataframe
#     return listOfPos
#
# print(wanted_table)
#
#
# '''To generate fasta of the scaffold containing the hits for further analysis'''
# # with open(tair_cds_fasta,'r') as tair, open(sinapis_fasta, 'r') as sinapis:
# #     outfile = open('/Users/josephbeegan/Documents/sinapis_scaffolds_with_matches.fa','w+')
# #     for sin_seq in SeqIO.parse(sinapis,'fasta'):
# #         sin_id = str(sin_seq.description)
# #         ids = wanted_table[wanted_table == sin_id].stack().index.tolist()
# #         for i in ids:
# #             wanted_id = wanted_table.loc[i[0],'sinapis_id']
# #             if wanted_id == sin_id:
# #                 print('WE DID IT')
# #                 outfile.write('>' + str(sin_id) + '\n'+ str(sin_seq.seq)+'\n')
#
# '''To generate file with the hits and sequences'''
# with open(tair_cds_fasta,'r') as tair, open(sinapis_fasta, 'r') as sinapis:
#     sin_seq_list = []
#     tair_seq_list = []
#     wanted_table["TAIR_SEQUENCE"] = np.nan
#     print(wanted_table)
#     for sin_seq in SeqIO.parse(sinapis,'fasta'):
#         sin_id = str(sin_seq.description)
#         ids = wanted_table[wanted_table == sin_id].stack().index.tolist()
#         if len(ids) != 0:
#             for i in ids:
#                 start_pos = wanted_table.loc[i[0],'sinapis_start']
#                 end_pos = wanted_table.loc[i[0],'sinapis_end']
#                 sin_seq_list.append(str(sin_seq.seq[start_pos:end_pos].upper()))
#
#     for tair_seq in SeqIO.parse(tair,'fasta'):
#         tair_id = str(tair_seq.description[0:9])
#         ids = wanted_table[wanted_table == tair_id].stack().index.tolist()
#         if len(ids) != 0:
#             for i in ids:
#                 tair_seq_list.append(wanted_table.loc[i[0],'tair_id'])
#                 wanted_table.loc[i[0],'TAIR_SEQUENCE'] = str(tair_seq.seq)
#
#
#     print(tair_seq_list)
#     wanted_table['SINAPIS_SEQUENCE'] = sin_seq_list
#
#
#
#
# print(wanted_table)
# wanted_table.to_csv('/Users/josephbeegan/Documents/sinapis_temp/default_blast_analysis/sinapis_to_tair_protein_analysis_results_dedup_with_sequences.txt', sep = '\t', index = False)
