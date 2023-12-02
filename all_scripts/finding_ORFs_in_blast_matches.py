import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd

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



reads = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/gene_prediction/RNAseq_to_nanopore_reads_10000_transcript.fa')
tair = fasta_condenser('/Volumes/sesame/joerecovery/genomes/TAIR10_cdna_20101214_updated')
cpc = fasta_condenser('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/gene_prediction/All_blast_to_TAIR10_cds_20101214_updated_dedup_parsed_AT2G46410_matches.fa')
print(cpc)

for i in range(len(cpc)):
    id = cpc['ID'][i]
    # print(id)
    # print(cpc['Seq'][i])
    if id[0] == 'A' and len(id) == 9:
        cpc['Seq'][i] = cpc['Seq'][i]
    else:
        location_index = reads.index[reads['ID']==id].to_list() #finding the index of the ID in the full reads file
        if len(location_index) > 1:
            print('error - duplicate sequence ids')
            break #Just a check to make sure the match is not empty
        ref_sequence = reads['Seq'][location_index[0]] #Extracting the full length read/sequence from the reads file
        loc_start = ref_sequence.find(cpc['Seq'][i]) 
        loc_end = loc_start + len(cpc['Seq'][i]) #Getting start and end points of blast hit in the reads file
        sequence = ref_sequence[loc_start -200:loc_end+200] #Adding bp either end of the blast hit to cast a wider net
        if len(sequence) % 3 != 0:
            sequence = sequence[0:len(sequence)-1]
            if len(sequence) % 3 != 0:
                sequence = sequence[0:len(sequence)-1]
        cpc['Seq'][i] = sequence #Replacing the blast hit sequence with the longer sequence from the extractio from the readfile



def find_orfs_with_trans(seq, trans_table, min_protein_length):
    print('RUNNING ORF FINDER NOW')
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        print(nuc)
        for frame in range(3):
            trans = nuc[frame:].translate(trans_table)
            print('TRANS:   ->' + trans)
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    answer.sort()
    print('ENDING ORF FINDER NOW')
    return answer


tair_orf_list = find_orfs_with_trans(Seq(cpc['Seq'][0]),11, 20)
sin_orf_list = find_orfs_with_trans(Seq(cpc['Seq'][1]),11, 20)

for start, end, strand, pro in tair_orf_list:
    print(
        "TAIR ORFS->>>>%s...%s - length %i, strand %i, %i:%i"
        % (pro[:30], pro[-3:], len(pro), strand, start, end)
    )
for start, end, strand, pro in sin_orf_list:
    print(
        "SINAPIS ORFS->>>>%s...%s - length %i, strand %i, %i:%i"
        % (pro[:30], pro[-3:], len(pro), strand, start, end)
    )