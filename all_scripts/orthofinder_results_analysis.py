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


def reverse_join_reversed_iter(s):
    s1 = ''.join(reversed(s))
    return s1
def extract_promoters(gff,genome,ID):
    trim_gff = gff[gff['ID'].str.contains(ID)]
    if trim_gff['strand'].iloc[0] == '+':
        trim_fasta_forward = genome[genome['ID'].isin(trim_gff['scaffold'])]
        print(trim_fasta_forward)
        trim_fasta_forward['set_seq'] = trim_fasta_forward['Seq'].str.slice(int(trim_gff['start'] -500),int(trim_gff['start']))
        trim_fasta_forward["seq_len"]= trim_fasta_forward["Seq"].str.len()
        trim_fasta_forward["set_seq_len"]= trim_fasta_forward["set_seq"].str.len()
        trim_fasta_forward['geneid'] = ID
        return trim_fasta_forward
    else:
        trim_fasta_reverse = genome[genome['ID'].isin(trim_gff['scaffold'])]
        print(trim_fasta_reverse)
        trim_fasta_reverse['set_seq'] = trim_fasta_reverse['Seq'].str.slice(int(trim_gff['end']),int(trim_gff['end'] + 500))
        trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(reverse_join_reversed_iter)
        trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(dna_complement)
        trim_fasta_reverse["seq_len"]= trim_fasta_reverse["Seq"].str.len()
        trim_fasta_reverse["set_seq_len"]= trim_fasta_reverse["set_seq"].str.len()
        trim_fasta_reverse['geneid'] = ID
        return trim_fasta_reverse
def dna_complement(x):
    seq = x.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    return seq
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
def dedup_list(list):
    seen=[]
    for i in range(len(list)):
        if list[i] not in seen:
            seen.append(list[i])
    return seen

sinapis_alba = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Phytozome/PhytozomeV13/Salba/v3.1/assembly/Salba_584_v3.0.fa')
arabidopsis = fasta_condenser('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_chr_all.fa')
sinapis_gff = gff_parser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.gene.gff3')
tair_gff = gff_parser('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff')
print(sinapis_gff)
print(tair_gff)


with open("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/orthofinder_hits.txt",'r') as f:
    tab = f.readlines()
    tab = [x.strip() for x in tab]


orthofinder_results=pd.read_csv("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/orthofinder_hits.txt")
print(orthofinder_results)
blast_results = pd.read_csv("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/BLASTs/Salba_584_v3.1.cds.fa_all_ids_BLAST_blast_to_TAIR10_cdna_20101214_updated_dedup.txt", sep='\t')
trim_blast_results = blast_results[['tair_id','sinapis_id']] 
pd.options.mode.chained_assignment = None


for n in range(len(orthofinder_results.index)):
    rown=orthofinder_results.iloc[n].to_list()
    print(rown[0])
    out_df=[]
    outname='/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/' + rown[0] + '_orthofinder_and_blast_hits_promoters_only.fa'   
    for i in range(len(rown)):
        try:
            if rown[i][0] != 'A':
                rown[i] = rown[i][0:15]
        except:
            continue
    sinapis_ids = blast_results['sinapis_id'][blast_results['tair_id'].str.contains(rown[0])].to_list()
    for j in range(len(sinapis_ids)):
        sinapis_ids[j] = sinapis_ids[j][0:15]
    sinapis_ids=dedup_list(sinapis_ids)
    full_list = dedup_list(rown +sinapis_ids)
    full_list = [x for x in full_list if str(x) != 'nan']
    for x in full_list:
        print(x)
        if x[0]=='A':
            df=extract_promoters(tair_gff,arabidopsis, x)
            print(('>' + df['geneid'].iloc[0] + '\n' + df['set_seq'].iloc[0] + '\n'))
        else:
            df=extract_promoters(sinapis_gff,sinapis_alba, x)
            print(('>' + df['geneid'].iloc[0] + '\n' + df['set_seq'].iloc[0] + '\n'))
        out_df.append(df)
    out_df = pd.concat(out_df)
    with open(outname,'w+') as sin_out:
        for i in range(len(out_df['ID'])):
            print(('>' + out_df['geneid'].iloc[i] + '\n' + out_df['set_seq'].iloc[i] + '\n'))
            sin_out.write('>' + out_df['geneid'].iloc[i] + '\n' + out_df['set_seq'].iloc[i] + '\n')






sinapis_alba_genes=fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/Salba_584_v3.0_genes_only_all.fa')
arabidopsis_genes=fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/tair_genes_only_all.fa')
all_genes_prom = pd.concat([sinapis_alba_genes, arabidopsis_genes], axis=0)

sinapis_alba_cds = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.cds.fa')
tair_cds = fasta_condenser('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_cds_20101214_updated')
all_genes = pd.concat([sinapis_alba_cds,tair_cds], axis=0)

 

# for i in range(len(tab)):
#     id_list = tab[i].split(',')
#     title=id_list[0]
#     arabidopsis_cds = tair_cds[tair_cds['ID'].str.contains(title)]
#     out_name = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/' + title + '_hits_and_orthofinder.fa'
#     with open(out_name,'w+') as out:
#         out.write('>' + arabidopsis_cds['ID'].iloc[0] + '\n' + arabidopsis_cds['Seq'].iloc[0] + '\n')
#         count = 0
#         seen = []
#         for id in id_list:
#             print(id[0:-2])
#             count += 1
#             found_gene = all_genes_prom[all_genes_prom['ID'].str.contains(id)]
#             if found_gene.empty:
#                 found_gene = all_genes_prom[all_genes_prom['ID'].str.contains(id[0:-2])]
#                 if found_gene.empty:
#                     found_gene = all_genes[all_genes['ID'].str.contains(id)]
#             if found_gene['ID'].iloc[0] not in seen:
#                 seen.append(found_gene['ID'].iloc[0])
#                 out.write('>' + found_gene['ID'].iloc[0] + '\n' + found_gene['Seq'].iloc[0] + '\n')
#             else:
#                 out.write('>' + found_gene['ID'].iloc[0] + "_" + str(count) + '\n' + found_gene['Seq'].iloc[0] + '\n')