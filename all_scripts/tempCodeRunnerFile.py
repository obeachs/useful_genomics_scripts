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

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def id_fix(id):
    if ',' in id:
        id = id.split(',')
        for i in range(len(id)):
            if 'Sialb' in id[i]:
                start = id[i].find('Sialb')
                id[i] = id[i][start:start+15]
            if 'Bra' in id[i]:
                start = id[i].find('Bra')
                id[i] = id[i][start:start +13]
            return ','.join(id)
    else:
        if 'Sialb' in id:
            start = id.find('Sialb')
            id = id[start:start+15]
        if 'Bra' in id:
            start = id.find('Bra')
            id = id[start:start +13]
        return id
def tair_id_fix(id):
    start = id.find('AT')
    return id[start:start+9]
def fix_alignment_num(string):
    start = string.find('-')
    newstring = string[:start]
    return int(newstring)
'''Setting up the analysis of the synteny between Sinapis alba (JGI v3) and Brassica rapa (v3)'''
'''Produces a dataframe of Syntenic blocks with the BLAST hits from the orthofinder output'''
'''
trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
trichome_ids_original = trichome_ids_original[0].to_list()
syntenic_blocks = pd.read_csv('/Users//josephbeegan/MCScanX/data/full_brapa_sinapis_edited.collinearity',sep='\t')
homologues = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/orthofinder/Results_Mar07/Orthogroups/Orthogroups.tsv', sep='\t')
homologues.rename(columns={'Salba_584_v3.1.protein':'Salba_584_v3'}, inplace=True)
homologues.rename(columns={'Brapa_genome_v3.0_pep':'Brapa_genome_v3'}, inplace=True)
homologues=homologues.assign(Araport11_pep_20220914=homologues['Araport11_pep_20220914'].str.split(',')).explode('Araport11_pep_20220914')
homologues=homologues.assign(Brapa_genome_v3=homologues['Brapa_genome_v3'].str.split(',')).explode('Brapa_genome_v3')
homologues=homologues.assign(Salba_584_v3=homologues['Salba_584_v3'].str.split(',')).explode('Salba_584_v3')
homologues=homologues[homologues['Araport11_pep_20220914'].notna()]
homologues['Araport11_pep_20220914'] = homologues['Araport11_pep_20220914'].apply(tair_id_fix)
homologues=homologues[homologues['Araport11_pep_20220914'].isin(trichome_ids_original)]
id_list = homologues['Brapa_genome_v3'].to_list() + homologues['Salba_584_v3'].to_list()
id_list = [x for x in id_list if str(x) != 'nan']
id_list = [id_fix(x) for x in id_list]

syntenic_blocks['Alignment'] = syntenic_blocks['Alignment'].apply(fix_alignment_num)
syntenic_groups = syntenic_blocks.groupby('Alignment')
count = syntenic_groups.ngroups
trichome_syntenic_blocks = []
for i in range(count):
    temp_df = syntenic_groups.get_group(i)
    if not temp_df[(temp_df['ID_1'].isin(id_list)) | (temp_df['ID_2'].isin(id_list))].empty:
        trichome_syntenic_blocks.append(temp_df)
syntenic_blocks = syntenic_blocks[(syntenic_blocks['ID_1'].isin(id_list)) | (syntenic_blocks['ID_2'].isin(id_list))]

trichome_syntenic_blocks = pd.concat(trichome_syntenic_blocks)
trichome_syntenic_blocks = trichome_syntenic_blocks[(trichome_syntenic_blocks['ID_1'].str.contains('Bra')) & (trichome_syntenic_blocks['ID_2'].str.contains('Sialb'))]
print(trichome_syntenic_blocks)
check_df = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/orthofinder/brapa_sinapis_alba_only/Results_Mar08/Orthogroups/Orthogroups.tsv',sep='\t')
check_df.rename(columns={'Salba_584_v3.1.protein':'Salba_584_v3'}, inplace=True)
check_df.rename(columns={'Brapa_genome_v3.0_pep':'Brapa_genome_v3'}, inplace=True)
check_df = check_df.dropna()
check_df['Brapa_genome_v3'] = check_df['Brapa_genome_v3'].apply(id_fix)
check_df['Salba_584_v3'] = check_df['Salba_584_v3'].apply(id_fix)
print(check_df)
blast_ID_1 = []
blast_ID_2 = []
for index, row in trichome_syntenic_blocks.iterrows():
    brapa_homologue= check_df[check_df['Salba_584_v3'].str.contains(row['ID_2'])]
    #trichome_syntenic_blocks.loc[index,'BLAST_HIT_ID1'] = brapa_homologue['Brapa_genome_v3'].to_list()
    blast_ID_2.append(brapa_homologue['Brapa_genome_v3'].to_list())
    sinapis_homologue= check_df[check_df['Brapa_genome_v3'].str.contains(row['ID_1'])]
    blast_ID_1.append(sinapis_homologue['Salba_584_v3'].to_list())
    #richome_syntenic_blocks.loc[index,'BLAST_HIT_ID2'] = sinapis_homologue['Salba_584_v3'].to_list()
trichome_syntenic_blocks['BLAST_HIT_ID2'] =blast_ID_2
trichome_syntenic_blocks['BLAST_HIT_ID1'] =blast_ID_1
trichome_syntenic_blocks = trichome_syntenic_blocks[['Alignment','ID_1','ID_2','BLAST_HIT_ID1','BLAST_HIT_ID2']]
print(len(blast_ID_2))
print(len(blast_ID_1))
#trichome_syntenic_blocks.to_csv('/Users/josephbeegan/MCScanX/data/trichome_syntenic_blocks_with_blast_hits.tsv', sep='\t', quoting=None,index=False)
'''
trichome_syntenic_blocks = pd.read_csv('/Users/josephbeegan/MCScanX/data/trichome_syntenic_blocks_with_blast_hits.tsv', sep='\t')
print(trichome_syntenic_blocks.columns)