import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
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

'''Make sure to have a '/' at the end of the directory string'''
directory = '/Volumes/seagatedrive/Project_folder/sinapis_assembly_shenanigans/blast_script_runs/SOAP_GSS_reads/'
full_sequence = '/Volumes/seagatedrive/Project_folder/sinapis_assembly_shenanigans/SOAP_SRR11653905_again.scafSeq'

summary_file_list = []
dedup_list = []

for f in listdir(directory):
    if 'dedup' in f:
        dedup_list.append(f)

if len(dedup_list) == 0:
    for f in listdir(directory):
        if 'short_summary' in f:
            summary_file_list.append(f)


def remove_dups(file,outfile):
    input_file = open(file,'r')
    lines_seen = set() # holds lines already seen
    with open(outfile, "w+") as output_file:
    	for each_line in input_file:
    	    if each_line not in lines_seen: # check if line is not duplicate
    	        output_file.write(each_line)
    	        lines_seen.add(each_line)




if len(summary_file_list) != 0:
    print('deduplicating files')
    for file in summary_file_list:
         file2 = directory + file
         newfilename = directory + file[:-4] + '_dedup.txt'
         remove_dups(file2,newfilename)

for f in listdir(directory):
    if 'dedup' in f:
        dedup_list.append(f)

id_tables = pd.read_table('/Volumes/seagatedrive/Project_folder/sinapis_assembly_shenanigans/tair_brapa_bnapus_bo_ids.txt',sep='\t')


''''Just showing what contigs matched to what reference genes'''

for i in range(len(dedup_list)):
    for f in dedup_list:
        if dedup_list[i] != f:
            with open(directory + dedup_list[i],'r') as dedup, open(directory + f,'r') as checkfile:
                for line in dedup:
                    for line2 in checkfile:
                        if line.split(':')[0] == line2.split(':')[0]:
                            print(line.split(':')[0] + ' matched to ' + line2.split(':')[1])

sinapis_ids_seen = []
for i in range(len(dedup_list)):
    for f in dedup_list:
        if dedup_list[i] != f:
            with open(directory + dedup_list[i],'r') as dedup, open(directory + f,'r') as checkfile:
                for line in dedup:
                    if line.split(' ')[0] not in sinapis_ids_seen:
                        print(line.split(' ')[0])
                        sinapis_ids_seen.append(line)

'''Create a unique list of the matches in the list'''
def get_all_matches_to_contig(contig,list):
    result_list = []
    for name in list:
        if contig in name:
            name_shortened = name.split(':')[1]
            name_shortened = name_shortened[:-1]
            result_list.append(name_shortened)
    seen = []
    for name in result_list:
        if name not in seen:
            seen.append(name)

    result_string = '_'.join(seen)
    return result_string


get_all_matches_to_contig('WIDR01028206',sinapis_ids_seen)



'''Taking the contigs that have more than 4 in the summary files, writing the title of the contigs and the
reference id that they matched and their sequence to a summary fasta file for analysis'''
summary_fasta_for_analysis = directory + 'contigs_with_several_hits.fasta'
with open(summary_fasta_for_analysis,'w+') as final_summary:
    for id in sinapis_ids_seen:
        count = 0
        for f in dedup_list:
            newid = id.split(':')[0]
            #newid = newid.split(' ')[0]
            newid = newid[:-3]

            with open (directory + f,'r') as file:
                for line in file:
                    if newid in line:
                        count += 1
        if count >= 4:
            print('DOOOPP DOOOP')
            for seq in SeqIO.parse(full_sequence,'fasta'):
                if seq.id in newid:
                    added_header = get_all_matches_to_contig(newid,sinapis_ids_seen)
                    # print('>' + seq.id + '_matched_to_' + added_header)
                    # print(seq.seq)
                    final_summary.write('>' + str(seq.id) + '_matched_to_' + added_header +'\n' + str(seq.seq) + '\n')
