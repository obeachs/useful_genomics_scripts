import Bio
import sys
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2
from Bio.Seq import Seq



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
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
def sinapis_id_fix(id):
    start = id.find('Sialb')
    return id[start:start+15]
def tair_id_fix(id):
    start = id.find('AT')
    return id[start:start+9]
hits = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/trichome_GRN_blast_hits_arabidopsis_singles_plus.tsv', sep='\t')
tair_prot = fasta_condenser('/Volumes/sesame//joerecovery//genomes/araport/Araport11_pep_20220914')
sinapis_prot = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.protein.fa')
brapa_prot = fasta_condenser('/Volumes//sesame/joerecovery/genomes/brapa/Brassica_rapa.Brapa_1.0.pep.all.fa')
trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
trichome_ids_original = trichome_ids_original[0].to_list()

seen =[]
count = 0
for index, row in hits.iterrows():
    count += 1
    fresh_df = []
    outitle = tair_id_fix(row['Araport11_pep_20220914'])
    outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/prot/raw_sequences/' + outitle  + '_trichome_GRN_homologues_proteins_with_brapa.fa'
    tair_df = tair_prot[tair_prot['ID'].str.contains(row['Araport11_pep_20220914'])]
    # print(tair_df)
    # print(row['Salba_584_v3'])
    # print(row['Brassica_rapa_pep'])0
    if row['Araport11_pep_20220914'] not in seen:
        seen.append(row['Araport11_pep_20220914'])
        fresh_df.append(tair_df)
        outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/prot/raw_sequences/' + outitle  + '_trichome_GRN_homologues_proteins_with_brapa.fa'
    else:
        outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/prot/raw_sequences/' + outitle  + '_' + str(count) + '_trichome_GRN_homologues_proteins_with_brapa.fa'
    fresh_df.append(tair_df)
    if isfloat(row['Salba_584_v3']) != True:
        sin_ids = row['Salba_584_v3'].split(',')
        for i in sin_ids:
            i = sinapis_id_fix(i)
            sin_df =sinapis_prot[sinapis_prot['ID'].str.contains(i)]
            fresh_df.append(sin_df)

    if isfloat(row['Brassica_rapa_pep']) != True:
        brapa_ids = row['Brassica_rapa_pep'].split(',')
        for i in brapa_ids:
            i = i[0:9]
            brapa_df =brapa_prot[brapa_prot['ID'].str.contains(i)]
            fresh_df.append(brapa_df)
    fresh_df = pd.concat(fresh_df)
    fresh_df = fresh_df.drop_duplicates()
    with open(outname,'w+') as sin_out:
        for i in range(len(fresh_df['ID'])):
            sin_out.write('>' + fresh_df['ID'].iloc[i] + '\n' + fresh_df['Seq'].iloc[i] + '\n')