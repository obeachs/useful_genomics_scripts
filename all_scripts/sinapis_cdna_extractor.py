import pandas as pd 
import Bio
from Bio import SeqIO
import numpy
import re
import subprocess
from glob import glob

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

def gff_parser(gff):
    readgff = pd.read_csv(gff, skiprows=3, sep='\t', header=None)
    readgff.columns = ['scaffold', 'annotation_origin','type','start','end','dot','strand','number','ID']
    readgff = readgff[readgff['type']=='gene']
    return readgff

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

def get_region(chr, start, end, ref):
    with open(ref,'r') as reference:
        for seq in SeqIO.parse(reference,'fasta'):
            if seq.id == chr:
                sequence = str(seq.seq[start:end])
    return sequence


salba_gff = gff_parser('/Volumes//sesame/joerecovery//Project_folder/sinapis_assembly_shenanigans//yang_assemblies//sinapis_alba/Sal.Chr.20210627.gff')
salba_gff['ID'] = salba_gff['ID'].apply(id_fix)
tair_gff = gff_parser('/Volumes//sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff')
tair_gff['ID'] = tair_gff['ID'].apply(id_fix)
salba = fasta_condenser('/Volumes//sesame/joerecovery//Project_folder/sinapis_assembly_shenanigans//yang_assemblies//sinapis_alba/Sal.Chr.20210627.fasta')
tair = fasta_condenser('/Volumes//sesame/joerecovery/genomes/TAIR/TAIR10_chr_all.fa')


print(get_region('Chr1',1, 10, '/Volumes//sesame/joerecovery/genomes/TAIR/TAIR10_chr_all.fa'))

#Incredibly dumb workaround to input the values from the gff dataframe into the script but it works


def write_cdna(fasta):
    out = fasta.split('.')[0] + '_full_length_genomic.fa'
    fasta = fasta_condenser(fasta)
    with open(out,'w+') as ou:
        for index, row in fasta.iterrows():
            if 'AT' in row['ID']:
                r = '/Volumes//sesame/joerecovery/genomes/TAIR/TAIR10_chr_all.fa'
                gff = tair_gff[tair_gff['ID'] == (row)['ID']]
                s = int(gff['start'])
                e = int(gff['end'])
                c = gff['scaffold'].iloc[0]
                print(c)
                seq = get_region(c,s,e,r)
                print(row['ID'])
                print(len(seq))
                ou.write('>' + str(row['ID']) + '\n' + seq + '\n')
            if 'Sal' in row['ID']:
                print(row['ID'])
                r = '/Volumes//sesame/joerecovery//Project_folder/sinapis_assembly_shenanigans//yang_assemblies//sinapis_alba/Sal.Chr.20210627.fasta'
                gff = salba_gff[salba_gff['ID'] == (row)['ID']]
                print(gff)
                s = int(gff['start'].iloc[0])
                e = int(gff['end'].iloc[0])
                c = gff['scaffold'].iloc[0]
                seq = get_region(c,s,e,r)
                print(str(row['ID']))
                print(len(seq))
                ou.write('>' + str(row['ID']) + '\n' + seq + '\n')


write_cdna('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/AT5G53200_pep_homologues.fa')
