import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir

cols = ['scaff','version','type','start','end','something','strand','full','ID']
def gff_parser(gff):
    readgff = pd.read_csv(gff, skiprows=3, sep='\t', header=None)
    readgff.columns = ['scaffold', 'annotation_origin','type','start','end','dot','strand','number','ID']
    readgff = readgff[readgff['type']=='gene']
    return readgff

def salba_id_fix(id):
    start = id.find('ID')
    newid = id[start+3:start+15]
    return newid
def tair_id_fix(id):
    start = id.find('AT')
    return id[start:start+9]
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

hits = pd.read_csv('/Users/josephbeegan/Salba_RNA/genelists/all_hits_strict_brapa_sep.csv')
print(hits)
brapa = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/promoters/brapa_promoters_only.fa')
print(brapa)
sinapis = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/promoters/Sal.Chr.20210627_promoters.fa')
tair = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/promoters/tair_promoters_only_all.fa')
tair['ID'] = tair['ID'].apply(tair_id_fix)

grouped = hits.groupby('tair')
for group_name, group_data in grouped:
    outfile = '/Users/josephbeegan/Salba_RNA/genelists/promoters' + group_name + '_promoter_homologs.fa'    
    group_data['gene_id'] = group_data['gene_id'].fillna('')
    group_data['tair'] = group_data['tair'].fillna('')
    group_data['brapa'] = group_data['brapa'].fillna('')
    print(group_name)
    print(group_data)
    # with open(outfile,'w+') as out:
    #     tair_check = tair[tair['ID'].str.contains(group_name)]
    #     if tair_check.empty:
    #         continue
    #     else:
    #         out.write('>' + tair_check['ID'].iloc[0] + '\n' + tair_check['Seq'].iloc[0] + '\n')
    #         for i in group_data['gene_id']:
    #             print(i)
    #             sin_check = sinapis[sinapis['ID'].str.contains(i)]
    #             print(sin_check['ID'])
    #             if sin_check.empty:
    #                 continue
    #             else:
    #                 out.write('>' + sin_check['ID'].iloc[0] + '\n' + sin_check['Seq'].iloc[0] + '\n')
    #         for i in group_data['brapa']:
    #             brapa_check = brapa[brapa['ID'].str.contains(i)]
    #             print(brapa_check['ID'])
    #             if brapa_check.empty:
    #                 continue
    #             else:
    #                 out.write('>' + brapa_check['ID'].iloc[0] + '\n' + brapa_check['Seq'].iloc[0] + '\n')