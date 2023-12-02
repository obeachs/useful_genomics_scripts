import pandas as pd
import numpy
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq


dict = {'A':'T','G':'C','T':'A','C':'G'}

test='ATCGATCGCTAGCTGACCCCC'
def reverse_join_reversed_iter(s):
    s1 = ''.join(reversed(s))
    return s1
print(reverse_join_reversed_iter(test))


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

alba_gff = gff_parser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.gene.gff3')
tair_gff = gff_parser('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff')
brapa_gff =  gff_parser('/Volumes/sesame/joerecovery/genomes/brapa/Brapa_genome_v3.0_genes.gff3')
print(brapa_gff)
with open("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt") as f:
    content_list = f.readlines()
    content_list = [x.strip() for x in content_list]

sinapis_alba = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Phytozome/PhytozomeV13/Salba/v3.1/assembly/Salba_584_v3.0.fa')
#inapis_flexuosa = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_arvensis/GWHBFWT00000000.genome.fasta')
arabidopsis = fasta_condenser('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_chr_all.fa')
brapa = fasta_condenser('/Volumes/sesame/joerecovery/genomes/brapa/Brapa_sequence_v3.0.fasta')
print(brapa)

pd.options.mode.chained_assignment = None
def extract_genomic_region(gff,genome,ID, promoters=0):
    if promoters==0:
        trim_gff = gff[gff['ID'].str.contains(ID)]
        trim_fasta = genome[genome['ID'].isin(trim_gff['scaffold'])]
        trim_fasta['set_seq'] = trim_fasta['Seq'].str.slice(int(trim_gff['start']) -200,int(trim_gff['end']) + 200)
        trim_fasta["seq_len"]= trim_fasta["Seq"].str.len()
        trim_fasta["set_seq_len"]= trim_fasta["set_seq"].str.len()
        trim_fasta['geneid'] = ID
        return trim_fasta
    else:
        trim_gff = gff[gff['ID'].str.contains(ID)]
        if trim_gff['strand'].iloc[0] == '+':
            trim_fasta_forward = genome[genome['ID'].isin(trim_gff['scaffold'])]
            trim_fasta_forward['set_seq'] = trim_fasta_forward['Seq'].str.slice(int(trim_gff['start'] -500),int(trim_gff['end']))
            trim_fasta_forward["seq_len"]= trim_fasta_forward["Seq"].str.len()
            trim_fasta_forward["set_seq_len"]= trim_fasta_forward["set_seq"].str.len()
            trim_fasta_forward['geneid'] = ID
            if trim_fasta_forward.empty:
                trim_gff = gff[gff['ID'].str.contains(ID)]
                trim_fasta_forward = genome[genome['ID'].isin(trim_gff['scaffold'])]
                trim_fasta_forward.loc['set_seq'] = trim_fasta_forward.loc['Seq'].str.slice(int(trim_gff['start']),int(trim_gff['end']))
                trim_fasta_forward["seq_len"]= trim_fasta_forward["Seq"].str.len()
                trim_fasta_forward["set_seq_len"]= trim_fasta_forward["set_seq"].str.len()
                trim_fasta_forward['geneid'] = ID
                return trim_fasta_forward
            return trim_fasta_forward
        else:
            trim_fasta_reverse = genome[genome['ID'].isin(trim_gff['scaffold'])]
            trim_fasta_reverse['set_seq'] = trim_fasta_reverse['Seq'].str.slice(int(trim_gff['start']),int(trim_gff['end'] + 500))
            trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(reverse_join_reversed_iter)
            trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(dna_complement)
            trim_fasta_reverse["seq_len"]= trim_fasta_reverse["Seq"].str.len()
            trim_fasta_reverse["set_seq_len"]= trim_fasta_reverse["set_seq"].str.len()
            trim_fasta_reverse['geneid'] = ID
            if trim_fasta_reverse.empty:
                trim_gff = gff[gff['ID'].str.contains(ID)]
                trim_fasta_reverse = genome[genome['ID'].isin(trim_gff['scaffold'])]
                trim_fasta_reverse.loc['set_seq'] = trim_fasta_reverse.loc['Seq'].str.slice(int(trim_gff['start']),int(trim_gff['end']))
                trim_fasta_reverse["seq_len"]= trim_fasta_reverse["Seq"].str.len()
                trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(reverse_join_reversed_iter)
                trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(dna_complement)
                trim_fasta_reverse["set_seq_len"]= trim_fasta_reverse["set_seq"].str.len()
                trim_fasta_reverse['geneid'] = ID
                return trim_fasta_reverse
            return trim_fasta_reverse

def extract_promoters(gff,genome,ID):
    trim_gff = gff[gff['ID'].str.contains(ID)]
    if trim_gff['strand'].iloc[0] == '+':
        trim_fasta_forward = genome[genome['ID'].isin(trim_gff['scaffold'])]
        if trim_fasta_forward.empty:
            trim_fasta_forward = genome[genome['ID'].str.contains(trim_gff['scaffold'].iloc[0])]
        trim_fasta_forward['set_seq'] = trim_fasta_forward['Seq'].str.slice(int(trim_gff['start'] -2500),int(trim_gff['start']))
        trim_fasta_forward["seq_len"]= trim_fasta_forward["Seq"].str.len()
        trim_fasta_forward["set_seq_len"]= trim_fasta_forward["set_seq"].str.len()
        trim_fasta_forward['geneid'] = ID
        return trim_fasta_forward
    else:
        trim_fasta_reverse = genome[genome['ID'].isin(trim_gff['scaffold'])]
        if trim_fasta_reverse.empty:
            trim_fasta_reverse = genome[genome['ID'].str.contains(trim_gff['scaffold'].iloc[0])]
        trim_fasta_reverse['set_seq'] = trim_fasta_reverse['Seq'].str.slice(int(trim_gff['end']),int(trim_gff['end'] + 2500))
        trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(reverse_join_reversed_iter)
        trim_fasta_reverse['set_seq'] = trim_fasta_reverse['set_seq'].apply(dna_complement)
        trim_fasta_reverse["seq_len"]= trim_fasta_reverse["Seq"].str.len()
        trim_fasta_reverse["set_seq_len"]= trim_fasta_reverse["set_seq"].str.len()
        trim_fasta_reverse['geneid'] = ID
        return trim_fasta_reverse



gff = gff_parser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.gene.gff3')
sinapis_matches = gff['ID'].str.slice(3,18).to_list()

tair_matches = tair_gff['ID'].str.slice(3,13).to_list()
brapa_mateches = brapa_gff['ID'].to_list()
# sinapis_genes_and_promoter = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/Salba_584_v3.0_genes_only_all.fa'
# tair_genes_and_promoter = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/tair_genes_only_all.fa'
#sinapis_promoter = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/Salba_584_v3.0_promoters_only.fa'
tair_promoter = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/promoters/tair_promoters_only_all.fa'
brapa_promoter= '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/promoters/brapa_promoters_only.fa'






sinapis_out_df = []
tair_out_df = []
brapa_out_df = []
# for i in sinapis_matches:
#     temp_df = extract_promoters(alba_gff, sinapis_alba,i)
#     sinapis_out_df.append(temp_df)
# sinapis_out_df = pd.concat(sinapis_out_df)

for i in tair_matches:
    # try:
    temp_df = extract_promoters(tair_gff,arabidopsis,i)
    # except:
    #     continue
    tair_out_df.append(temp_df)
tair_out_df = pd.concat(tair_out_df)
print(tair_out_df)

# for i in brapa_mateches:
#     temp_df = extract_promoters(brapa_gff, brapa,i)
#     brapa_out_df.append(temp_df)
# brapa_out_df = pd.concat(brapa_out_df)
# print(brapa_out_df)

# with open(sinapis_promoter,'w+') as sin_out:
#     for i in range(len(sinapis_out_df['ID'])):
#         sin_out.write('>' + sinapis_out_df['geneid'].iloc[i] + '\n' + sinapis_out_df['set_seq'].iloc[i] + '\n')

# with open(brapa_promoter,'w+') as brapa_out:
#     for i in range(len(brapa_out_df['ID'])):
#         brapa_out.write('>' + brapa_out_df['geneid'].iloc[i] + '\n' + brapa_out_df['set_seq'].iloc[i] + '\n')

with open(tair_promoter,'w+') as tair_out:
    for i in range(len(tair_out_df['ID'])):
        tair_out.write('>' + tair_out_df['geneid'].iloc[i] + '\n' + tair_out_df['set_seq'].iloc[i] + '\n')
