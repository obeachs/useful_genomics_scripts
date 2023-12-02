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
parser.add_argument('--refcdna', required = False, default = '/Volumes/sesame/joerecovery/genomes/araport/Araport11_pep_20220914', type = str)
parser.add_argument('--queryfasta',required = False, default = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Phytozome/PhytozomeV13/Salba/v3.1/assembly/Salba_584_v3.0.fa', type = str)
parser.add_argument('--blastoutput', required = True, default = None, type= str)
parser.add_argument('--geneids', required = False, default = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', type= str)
parser.add_argument('--noextend', required=False, default=0,type=int)
parser.add_argument('--allgenes', required=False, default=0, type=int)
args = parser.parse_args()

#This just works a lot faster than SeqIO.parse - especially because we don't know the number
#of sequences
#Adding tair = 1 will just change the titles of the columns that are created in the final 
#table
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




#Reading in the deduplicated blast analysis output and removing characters
#after the decimal point - this is mostly just to get hits to genes as 
#opposed to different transcripts
table = pd.read_table(args.blastoutput)
table['tair_id'] =table['tair_id'].replace(r'\..*', '', regex=True)
#Reading in the fastas of the blast reference and query
tair_fasta = args.refcdna
sinapis_genome = args.queryfasta
geneids = args.geneids
outfile = args.blastoutput
if args.allgenes != 0:
    outfile = outfile.split('.')[0] + '_parsed_all_genes.fa'
outfile = outfile.split('.')[0] + '_parsed.fa'
sinapis = fasta_condenser(sinapis_genome)
tair = fasta_condenser(tair_fasta,tair=1)
print(tair)


if '/' in geneids:
    with open(geneids,'r') as ids:
        id_list = ids.read().splitlines() 
else:
    id_list = []
    geneids = geneids
    id_list.append(geneids)
#Condensing the TAIR/reference fasta to just be genes we are interested in
if args.allgenes == 0:
    tri_tair = tair[tair['TAIR_ID'].isin(id_list)]
    matches_table = table[table['tair_id'].isin(id_list)]
    print(matches_table)
    #Condensing the query fasta to just be genes we are interested in
    tri_sinapis = sinapis[sinapis['ID'].isin(matches_table['sinapis_id'])]
    #Merging the tables so that the full sequences of the hits can be extracted
    merged = matches_table.merge(tri_tair, how='left',left_on='tair_id', right_on='TAIR_ID')
    merged = merged.merge(tri_sinapis, how='left',left_on='sinapis_id', right_on='ID')
    #Unsure why this next line is needed - sometimes new query IDs seem to be added
    # - will try to fix
    merged = merged[merged['Seq'].notna()]

    #TAIR_Seq_x is the cdna, y is the cds
    print(merged)
    #print(merged['Seq'])
    #I am sure there is a faster way to do this but just converting each column
    #to a list to make it easier to loop through them
    sinapis_id = merged['sinapis_id'].to_list()
    sinapis_seq = merged['Seq'].to_list()
    print(sinapis_seq)
    sinapis_start = merged['sinapis_start'].to_list()
    sinapis_end = merged['sinapis_end'].to_list()
    tair_ids = merged['tair_id'].to_list()
    tair_seq_cdna = merged['TAIR_Seq'].to_list()
else:
    matches_table = table
    print(matches_table)
    #Condensing the query fasta to just be genes we are interested in
    #tri_tair = tair['TAIR_ID']
    #tri_sinapis = sinapis[sinapis['ID'].isin(matches_table['sinapis_id'])]
    tri_sinapis = sinapis
    tri_tair = tair
    #Merging the tables so that the full sequences of the hits can be extracted
    merged = matches_table.merge(tri_tair, how='left',left_on='tair_id', right_on='TAIR_ID')
    merged = merged.merge(tri_sinapis, how='left',left_on='sinapis_id', right_on='ID')
    merged['set_seq'] = merged['Seq'].str.slice(merged['sinapis_start'],merged['sinapis_end'])
    #Unsure why this next line is needed - sometimes new query IDs seem to be added
    # - will try to fix
    merged = merged[merged['Seq'].notna()]
    merged['sinapis_start'] = merged['sinapis_start']-1
    print(merged.columns)
    print(merged['TAIR_Seq'])
    print(merged['Seq'])
    #TAIR_Seq_x is the cdna, y is the cds
    #print(merged['Seq'])
    #I am sure there is a faster way to do this but just converting each column
    #to a list to make it easier to loop through them
    sinapis_id = merged['sinapis_id'].to_list()
    sinapis_seq = merged['Seq'].to_list()
    sinapis_start = merged['sinapis_start'].to_list()
    sinapis_end = merged['sinapis_end'].to_list()
    tair_ids = merged['tair_id'].to_list()
    tair_seq_cdna = merged['TAIR_Seq'].to_list()
    
#No extend option should be used for transcripts mostly, just because otherwise the 
#index will be way out of bounds
'''
for i in range(len(sinapis_seq)):
    try:
        sinapis_seq[i] = sinapis_seq[i][sinapis_start[i]-500:sinapis_end[i]-500]
    except:
        sinapis_seq[i] = sinapis_seq[i]
'''

with open(outfile,'w+') as out:
    for i in range(len(tair_ids)):
        out.write('>' + tair_ids[i]+'_cdna' + '\n' + tair_seq_cdna[i] + '\n' +
         '>' + tair_ids[i] +'_cds' + '\n' + '>' + sinapis_id[i] + '\n' + 
         sinapis_seq[i] + '\n')





