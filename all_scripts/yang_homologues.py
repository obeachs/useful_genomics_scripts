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
from Bio.Seq import Seq
from os import listdir

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
def tair_id_fix(id):
    start = id.find('AT')
    return id[start:start+9]
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
def id_fix(id):
    if 'Sialb' in id:
        start = id.find('Sialb')
        id = id[start:start+15]
    if 'Bra' in id:
        start = id.find('Bra')
        id = id[start:start +13]
    if 'AT' in id:
        start = id.find('AT')
        id = id[start:start + 9]
    if 'Sal' in id:
        start = id.find('Sal')
        id = id[start:start+12]
    return id


blocks = pd.read_csv('/Volumes/sesame//joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/orthofinder_results_yang_tair_brapa/Orthogroups/Orthogroups.tsv',sep='\t')
blocks = blocks[blocks['Araport11_pep_20220914'].notna()]
blocks['Araport11_pep_20220914'] = blocks['Araport11_pep_20220914'].apply(tair_id_fix)
sin_cds = fasta_condenser('/Volumes//sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Sal.cds')
sin_pep = fasta_condenser('/Volumes//sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Sal.pep')
brapa_cds = fasta_condenser('/Volumes/sesame/joerecovery//genomes//brapa/Brapa_genome_v3.0_cds.fasta')
brapa_pep = fasta_condenser('/Volumes/sesame/joerecovery//genomes//brapa/Brapa_genome_v3.0_pep.fasta')

tair_cds = fasta_condenser('/Volumes//sesame//joerecovery//genomes/TAIR/TAIR10_cds_20101214_updated')
tair_pep = fasta_condenser('/Volumes/sesame/joerecovery/genomes/araport/Araport11_pep_20220914')
trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
trichome_ids_original = trichome_ids_original[0].to_list()

blocks = blocks[blocks['Araport11_pep_20220914'].isin(trichome_ids_original)]
print(blocks)
# for index, row in blocks.iterrows():
#     tair_id = row['Araport11_pep_20220914']
#     tair_seq = tair_cds['Seq'][tair_cds['ID'].str.contains(tair_id)].to_list()[0]
#     sin_ids = []
#     sin_seqs = []
#     if isfloat(row['Salba_yang_pep']) != True:
#         for i in row['Salba_yang_pep'].split(','):
#             i = id_fix(i)
#             sin_ids.append(i)
#             hom = sin_cds[sin_cds['ID'].str.contains(i)]
#             seq = hom['Seq'].to_list()[0]
#             sin_seqs.append(seq)
#     brapa_ids = []
#     brapa_seqs = []
#     if isfloat(row['Brapa_genome_v3.0_pep']) != True:
#         for i in row['Brapa_genome_v3.0_pep'].split(','):
#             i = id_fix(i)
#             brapa_ids.append(i)
#             hom = brapa_cds[brapa_cds['ID'].str.contains(i)]
#             seq = hom['Seq'].to_list()[0]
#             brapa_seqs.append(seq)
#     print(brapa_ids)
#     print(brapa_seqs)
#     print(sin_ids)
#     print(sin_seqs)
#     outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/' + tair_id + '_cds_homologues.fa'
#     with open(outname,'w+') as out:
#         out.write('>' + str(tair_id) + '\n' + str(tair_seq) + '\n')
#         for sid,seq in zip(sin_ids, sin_seqs):
#             out.write('>' + str(sid) + '\n' + str(seq) + '\n')
#         for brid,beq in zip(brapa_ids, brapa_seqs):
#             out.write('>' + str(brid) + '\n' + str(beq) + '\n')

for index, row in blocks.iterrows():
    tair_id = row['Araport11_pep_20220914']
    tair_seq = tair_pep['Seq'][tair_pep['ID'].str.contains(tair_id)].to_list()[0]
    sin_ids = []
    sin_seqs = []
    if isfloat(row['Salba_yang_pep']) != True:
        for i in row['Salba_yang_pep'].split(','):
            i = id_fix(i)
            sin_ids.append(i)
            hom = sin_pep[sin_pep['ID'].str.contains(i)]
            seq = hom['Seq'].to_list()[0]
            sin_seqs.append(seq)
    brapa_ids = []
    brapa_seqs = []
    if isfloat(row['Brapa_genome_v3.0_pep']) != True:
        for i in row['Brapa_genome_v3.0_pep'].split(','):
            i = id_fix(i)
            brapa_ids.append(i)
            hom = brapa_pep[brapa_pep['ID'].str.contains(i)]
            seq = hom['Seq'].to_list()[0]
            brapa_seqs.append(seq)
    print(brapa_ids)
    print(brapa_seqs)
    print(sin_ids)
    print(sin_seqs)
    outname = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/' + tair_id + '_pep_homologues.fa'
    with open(outname,'w+') as out:
        out.write('>' + str(tair_id) + '\n' + str(tair_seq) + '\n')
        for sid,seq in zip(sin_ids, sin_seqs):
            out.write('>' + str(sid) + '\n' + str(seq) + '\n')
        for brid,beq in zip(brapa_ids, brapa_seqs):
            out.write('>' + str(brid) + '\n' + str(beq) + '\n')
