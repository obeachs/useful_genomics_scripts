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

trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
trichome_ids_original = trichome_ids_original[0].to_list()
sinapis_pep = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.protein.fa')
brapa_pep = fasta_condenser('/Volumes/sesame/joerecovery/genomes/brapa/Brapa_genome_v3.0_pep.fasta')
orthogroups = pd.read_csv('/Volumes//sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/orthofinder/Results_Mar07/Orthogroups/Orthogroups.tsv',sep='\t')
orthogroups.rename(columns={'Salba_584_v3.1.protein':'Salba_584_v3'}, inplace=True)
orthogroups.rename(columns={'Brapa_genome_v3.0_pep':'Brapa_genome_v3'}, inplace=True)
orthogroups=orthogroups.assign(Brapa_genome_v3=orthogroups['Brapa_genome_v3'].str.split(',')).explode('Brapa_genome_v3')
orthogroups=orthogroups.assign(Salba_584_v3=orthogroups['Salba_584_v3'].str.split(',')).explode('Salba_584_v3')
orthogroups=orthogroups.assign(Araport11_pep_20220914=orthogroups['Araport11_pep_20220914'].str.split(',')).explode('Araport11_pep_20220914')
orthogroups = orthogroups[orthogroups['Araport11_pep_20220914'].notna()]
orthogroups['Araport11_pep_20220914'] = orthogroups['Araport11_pep_20220914'].apply(tair_id_fix)
orthogroups = orthogroups[orthogroups['Araport11_pep_20220914'].isin(trichome_ids_original)]
print(orthogroups['Salba_584_v3'])
sin_outfasta = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/salba_trichome_protein.fasta'
brapa_outfasta = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/brapa_trichome_protein.fasta'
with open(sin_outfasta, 'w+') as sin_fasta, open(brapa_outfasta,'w+') as brapa_fasta:
    for index, row in orthogroups.iterrows():
        if isfloat(row['Salba_584_v3']) != True:
            sin_id = row['Salba_584_v3']
            sin_out = sinapis_pep[sinapis_pep['ID'].str.contains(sin_id)]
            if not sin_out.empty:
                sin_fasta.write('>' + sin_id + '\n' + sin_out['Seq'].iloc[0] +'\n')
        if isfloat(row['Brapa_genome_v3']) != True:
            brapa_id = row['Brapa_genome_v3']
            brapa_out = brapa_pep[brapa_pep['ID'].str.contains(brapa_id)]
            if not brapa_out.empty:
                brapa_fasta.write('>' + brapa_id + '\n' + brapa_out['Seq'].iloc[0] +'\n')