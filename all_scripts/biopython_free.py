import Bio
import sys
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio.Seq import Seq

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
 
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

def gff_parser(gff):
    readgff = pd.read_csv(gff, skiprows=3, sep='\t', header=None)
    readgff.columns = ['scaffold', 'annotation_origin','type','start','end','dot','strand','number','ID']
    readgff = readgff[readgff['type']=='gene']
    return readgff

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

# cdsfile = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.cds.fa')
# id = 'Sialb.0255s0110.1.p'
# out = cdsfile['Seq'][cdsfile['ID'].str.contains(id)]
# print(out)


trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
trichome_ids_original = trichome_ids_original[0].to_list()
ortho = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Results_Mar14/Orthogroups/Orthogroups.tsv',sep='\t')
ortho = ortho[ortho['Araport11_pep_20220914'].notna()]
ortho =ortho.assign(Araport11_pep_20220914=ortho['Araport11_pep_20220914'].str.split(',')).explode('Araport11_pep_20220914')
ortho['Araport11_pep_20220914'] = ortho['Araport11_pep_20220914'].apply(id_fix)
out = ortho[ortho['Araport11_pep_20220914'].isin(trichome_ids_original)]
print(ortho['Araport11_pep_20220914'])
print(trichome_ids_original)
out.to_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Results_Mar14/Orthogroups/Orthogroups_trichomes.tsv',sep='\t',quoting=None,index=False)

# wd = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Results_Mar14/WorkingDirectory/'
# blast_files_yang = [wd + 'Blast1_1.txt',wd+'Blast1_2.txt', wd +'Blast2_1.txt',wd+'Blast2_2.txt']
# blast_files = list_full_paths(wd)



# brapa_fa = '/Volumes/sesame/joerecovery/genomes/brapa/Brapa_genome_v3.0_pep.fasta'
# sinapis_fa = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Sal.pep'
# id_swaps = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Results_Mar14/WorkingDirectory/SequenceIDs.txt',sep='\t', header=None)
# print(id_swaps)
# id_swaps.columns = ['Ortho','ID']
# id_swaps['ID'] = id_swaps['ID'].apply(id_fix)
# id_swaps_dict = id_swaps.set_index('Ortho')['ID'].to_dict()

# for i in blast_files_yang:
#     if 'Blast' in i:
#         outname = i.split('.')[0] + '_edited.txt'
#         ortho = pd.read_csv(i, sep='\t', header=None)
#         print(ortho)
#         ortho.columns = ['name_1', 'name_2', 'score', 'num_1','num_2','num_3','num_4','num_5','num_6','num_7','num_8', 'num_9']
#         ortho['name_1'] = ortho['name_1'].map(id_swaps_dict)
#         ortho['name_2'] = ortho['name_2'].map(id_swaps_dict)
#         ortho.to_csv(outname, quoting=None,index=False,header=None)

# sinapis_gff = gff_parser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Sal.Chr.20210627.gff')
# brapa_gff = gff_parser('/Volumes/sesame/joerecovery/genomes/brapa/Brapa_genome_v3.0_genes.gff3')
# sinapis_gff = sinapis_gff[['scaffold','ID','start','end']]
# sinapis_gff['scaffold'] = sinapis_gff['scaffold'].str.replace('scaffold','sinapis_scaff')
# sinapis_gff['ID'] = sinapis_gff['ID'].apply(id_fix)
# brapa_gff = brapa_gff[['scaffold','ID','start','end']]
# brapa_gff['ID'] = brapa_gff['ID'].apply(id_fix)

# brapa_gff.to_csv('/Users/josephbeegan/MCScanX/data/brapa.gff', quoting=None,header=None,index=False,sep='\t')
# sinapis_gff.to_csv('/Users/josephbeegan/MCScanX/data/sinapis_new.gff', quoting=None,header=None,index=False, sep='\t')






# def tair_id_fix(id):
#     start = id.find('AT')
#     return id[start:start+9]
# grn = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/trichome_GRN_blast_hits_arabidopsis_singles.tsv', sep='\t')
# list = grn['Araport11_pep_20220914'].to_list()
# print(len(list))
# trichome_ids_original = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', header=None)
# trichome_ids_original = trichome_ids_original[0].to_list()
# list = list + trichome_ids_original
# seen = []
# for i in list:
#     i = tair_id_fix(i)
#     if i not in seen:
#         seen.append(i)
#     else:
#         continue
# list = seen
# print(len(list))
# tab = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/orthofinder/Results_Mar03/Orthogroups/Orthogroups.tsv',sep='\t')
# tab_tair_only =tab.assign(Araport11_pep_20220914=tab['Araport11_pep_20220914'].str.split(',')).explode('Araport11_pep_20220914')
# tab_tair_only = tab_tair_only[tab_tair_only['Araport11_pep_20220914'].notna()]
# tab_tair_only.rename(columns={'Salba_584_v3.1.protein':'Salba_584_v3'}, inplace=True)
# tab_tair_only.rename(columns={'Brassica_rapa.Brapa_1.0.pep.all':'Brassica_rapa_pep'}, inplace=True)
# tab_tair_only['Araport11_pep_20220914'] = tab_tair_only['Araport11_pep_20220914'].apply(tair_id_fix)
# tab_tair_only = tab_tair_only[tab_tair_only['Araport11_pep_20220914'].isin(list)]
# print(len(tab_tair_only['Araport11_pep_20220914']))
#tab_tair_only.to_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/GRN/trichome_GRN_blast_hits_arabidopsis_singles_plus.tsv', sep='\t')

# genelist = []
# with open('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/trichome_GRN_list.txt','r') as out:
#     for line in out:
#         genelist.append(line.strip())

# df = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/trichome_GRN_blast_hits.tsv', sep='\t')
# df.rename(columns={'Salba_584_v3.1.protein':'Salba_584_v3'}, inplace=True)
# df = df.assign(Araport11_pep_20220914=df['Araport11_pep_20220914'].str.split(',')).explode('Araport11_pep_20220914')
# df['Araport11_pep_20220914'].str.strip()
# df.to_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/trichome_GRN_blast_hits_arabidopsis_singles.tsv', sep='\t', index=False)



# def listToString(s):
 
#     # initialize an empty string
#     str1 = ""
 
#     # traverse in the string
#     for ele in s:
#         str1 += ele
 
#     # return string

#     return str1
# zainab_ortho = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/zainab/Orthogroups.tsv', sep = '\t')
# sin_gtf = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/zainab/SRR2961890_Salba_new.gtf', sep = '\t',header=None)
# sin_gtf[8] = sin_gtf[8].str.replace('gene_id ','')
# sin_gtf[8] = sin_gtf[8].str.replace('"','')
# sin_gtf[9] = sin_gtf[9].str.replace('transcript_id ','')
# sin_gtf[9] = sin_gtf[9].str.replace('"','')
# ortholist = zainab_ortho['SRR2961890_Salba_new_augustus_prot'].to_list()
# print(ortholist[12197])
# if 'g39560.t1' in ortholist[12197]:
#     print('what')
# gtf_list = sin_gtf[9].to_list()
# ortholist_string =[]
# for i in gtf_list:
#     try:
#         print(zainab_ortho[zainab_ortho['SRR2961890_Salba_new_augustus_prot'].str.contains(i)])
#     except:
#         continue

# print(ortholist_string)

# for i in range(len(gtf_list)):
#     print(gtf_list[i])
#     if gtf_list[i] in ortholist:
#         print('FARST')



# with open('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.protein.fa', 'r') as genome:
#     count=0
#     for seq in SeqIO.parse(genome,'fasta'):
#         if '.1.p' in seq.description:
#             count += 1

#     print(count)

# def fasta_condenser(fasta, tair=0):
#     '''Preferably use this on the fasta file that has less info/annotation
#     with it - usually the query fasta'''
#     namelist = []
#     seqlist = []
#     if tair==0:
#         with open (fasta,'r') as fa:
#             for seq in SeqIO.parse(fa,'fasta'):
#                 namelist.append(seq.description)
#                 seqlist.append(str(seq.seq))
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
# # print(reads)
#     else:
#         with open (fasta,'r') as fa:
#             for seq in SeqIO.parse(fa,'fasta'):
#                 namelist.append(seq.id[0:9])
#                 seqlist.append(str(seq.seq))
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
#     return df




# # reads = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_flexuosa_nanopore/all_runs_to_jan_2023/all_runs_combined.fa')


# # reads['length']  = reads['Seq'].str.len()
# # print(reads)
# # print(reads.sort_values('length', ascending=False))
