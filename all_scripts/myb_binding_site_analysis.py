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
import glob
from glob import glob
from os import listdir
from Bio import pairwise2
from Bio.Seq import Seq

def dna_complement(x):
    seq = x.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    seq = seq[::-1]
    return seq

def Convert(string):
    list1 = []
    list1[:0] = string
    return list1
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
def add_regex_column_check(expression, df, colname):
    df['Reverse_comp_seq'] = df['Seq'].apply(dna_complement)
    new_column_forward = []    
    for i in range(len(df['ID'])):
        if len(re.findall(expression, df['Seq'].iloc[i])) !=0:
            new_column_forward.append('Yes')
        elif len(re.findall(expression, df['Reverse_comp_seq'].iloc[i])) !=0:
            new_column_forward.append('Yes')
        else:
            new_column_forward.append('No')
    df[colname] = new_column_forward

trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
trichome_ids_original = trichome_ids_original[0].to_list()

myb_core = r'[CT][AGCT]GTT[GT]'
mbs_1 = r'C[ACGT]GTT[AG]'
mbs_2 = r'G[GT]T[AT]GTT[AG]'
mbs_3 = r'[CT]ACC[AT]A[AC]C'
r2_myb_1 = r'[CT]AAC[AGCT]G'
r2_myb_2 = r'ACC[AT]A[AC]'

def sinapis_id_fix(id):
    start = id.find('Sialb')
    return id[start:start+15]
def tair_id_fix(id):
    start = id.find('AT')
    return id[start:start+9]

genelist_df = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/trichome_GRN_blast_hits_arabidopsis_singles_plus.tsv',sep='\t')

genelist_df['Araport11_pep_20220914'] =genelist_df['Araport11_pep_20220914'].apply(tair_id_fix)
tair_genelist = genelist_df['Araport11_pep_20220914'].to_list()


tair_promoters='/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/tair_promoters_only_all.fa'
promoters = fasta_condenser(tair_promoters)
promoters['ID'] = promoters['ID'].apply(tair_id_fix)
promoters = promoters[promoters['ID'].isin(tair_genelist)]

test_df = genelist_df[genelist_df['Araport11_pep_20220914'].isin(trichome_ids_original)]
print(test_df)

sinapis_promoters = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/Salba_584_v3.0_promoters_only.fa'
sin_promoters = fasta_condenser(sinapis_promoters)
temp_df_for_sinapis_id_expansion =genelist_df.assign(Salba_584_v3=genelist_df['Salba_584_v3'].str.split(',')).explode('Salba_584_v3')
temp_df_for_sinapis_id_expansion = temp_df_for_sinapis_id_expansion[temp_df_for_sinapis_id_expansion['Salba_584_v3'].notna()]
sinapis_genelist = temp_df_for_sinapis_id_expansion['Salba_584_v3'].to_list()
sinapis_genelist = [sinapis_id_fix(elem) for elem in sinapis_genelist]
print(sinapis_genelist)

brapa_promoter = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/brapa_promoters_only_all.fa')
temp_df_for_brapa_id_expansion =genelist_df.assign(Brassica_rapa_pep=genelist_df['Brassica_rapa_pep'].str.split(',')).explode('Brassica_rapa_pep')
temp_df_for_brapa_id_expansion = temp_df_for_brapa_id_expansion[temp_df_for_brapa_id_expansion['Brassica_rapa_pep'].notna()]
brapa_genelist = temp_df_for_brapa_id_expansion['Brassica_rapa_pep'].to_list()
brapa_genelist = [elem[:9] for elem in brapa_genelist]

'''Final dataframes containing all of the genes/homologues involved in trichome GRN'''
tair_promoters_sub =promoters[promoters['ID'].isin(tair_genelist)]
print(sin_promoters['ID'])
sinapis_promoters_sub = sin_promoters[sin_promoters['ID'].isin(sinapis_genelist)]
print(sinapis_promoters_sub)
brapa_promoters_sub = brapa_promoter[brapa_promoter['ID'].isin(brapa_genelist)]
print(brapa_promoters_sub)

add_regex_column_check(myb_core, promoters,'myb_core')
add_regex_column_check(mbs_1, promoters,'mbs_1')
add_regex_column_check(mbs_2, promoters,'mbs_2')
add_regex_column_check(mbs_2, promoters,'mbs_3')
'''This section is for the entire trichome GRN in arbaidopsis'''


'''This section is for the major players in the MBW trichome complex'''
trichome_promoters = []
for i in glob("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/*fa"):
    title = i.split('.')[0]
    title = title+ 'MYB_site_check.csv'
    trichome_promoters.append(fasta_condenser(i))
for i in trichome_promoters:
    add_regex_column_check(myb_core, i,'myb_core')
    add_regex_column_check(mbs_1, i,'mbs_1')
    add_regex_column_check(mbs_2, i,'mbs_2')
    add_regex_column_check(mbs_3, i,'mbs_3')
    add_regex_column_check(r2_myb_1, i,'r2_myb_1')
    add_regex_column_check(r2_myb_2, i,'r2_myb_2')

for i in trichome_promoters:
    if i['ID'].iloc[0] in trichome_ids_original:
        print(i['ID'].iloc[0])
    title = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/TFBS/' + i['ID'].iloc[0] + '_MYB_BS_check_brapa.csv'
    i.to_csv(title)
'''THIS WORKS'''
'''
seen =[]
count = 0
for index, row in genelist_df.iterrows():
    count += 1
    fresh_df = []
    tair_df = tair_promoters_sub[tair_promoters_sub['ID'] == row['Araport11_pep_20220914']]
    if row['Araport11_pep_20220914'] not in seen:
        seen.append(row['Araport11_pep_20220914'])
        fresh_df.append(tair_df)
        outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/' + row['Araport11_pep_20220914'] + '_trichome_GRN_homologues_promoters_with_brapa.fa'
    else:
        tair_df['ID'] = tair_df['ID'] + '_' + str(count)
        fresh_df.append(tair_df)
        outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/' + row['Araport11_pep_20220914'] + '-' + str(count) + '_trichome_GRN_homologues_promoters_with_brapa.fa'
    if isfloat(row['Salba_584_v3']) != True:
        sin_ids = row['Salba_584_v3'].split(',')
        for i in sin_ids:
            i = sinapis_id_fix(i)
            sin_df = sinapis_promoters_sub[sinapis_promoters_sub['ID']==i]
            fresh_df.append(sin_df)

    if isfloat(row['Brassica_rapa_pep']) != True:
        brapa_ids = row['Brassica_rapa_pep'].split(',')
        for i in brapa_ids:
            i = i[0:9]
            brapa_df = brapa_promoters_sub[brapa_promoters_sub['ID']==i]
            fresh_df.append(brapa_df)
    fresh_df = pd.concat(fresh_df)
    fresh_df = fresh_df.drop_duplicates()
    with open(outname,'w+') as sin_out:
        for i in range(len(fresh_df['ID'])):
            sin_out.write('>' + fresh_df['ID'].iloc[i] + '\n' + fresh_df['Seq'].iloc[i] + '\n')
'''