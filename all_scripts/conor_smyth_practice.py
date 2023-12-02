import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
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

sequence = 'GAACTCTCTGTTCTTTGTGATGCTGAAGTCACACTCATCGTCTTCTCTAGACCGTTTGATTCACGGAACTATTTCCAAGTCGCTGCATTGCAACCTAACAATCACCATTACTCATCCGCAGGTCGCCAAGACCAAACCGCTCTTCAGTTAGTGTAA'

gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def translate_dna(dna):
    last_codon_start = len(dna) - 2
    protein = ""
    for start in range(0,last_codon_start,3):
        codon = dna[start:start+3]
        aa = gencode.get(codon.upper(), 'X')
        protein = protein + aa
    return protein


print(translate_dna(sequence))

conor_tair_to_bnapus = '/Users/josephbeegan/Desktop/Bnapus_to_TAIR_6targets_xml'

id_tables = pd.read_table('/Volumes/seagatedrive/sinapis_assembly_shenanigans/tair_brapa_bnapus_bo_ids.txt',sep='\t')
tair_ids = id_tables['TAIR'].to_list()


'''general for loop structure for BLAST analysis'''
with open(conor_tair_to_bnapus,'r') as conor_tair_to_bnapus_open:
    parsed_file = NCBIXML.parse(conor_tair_to_bnapus_open)
    for record in parsed_file:
        for description in record.descriptions:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    '''Conditional arguments and writing to files goes here'''




print('hi')

        print(line)



biopythonsequence = Seq(sequence)
print(biopythonsequence)
print(biopythonsequence.reverse_complement())


with open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/new_blast_fasta_brapa_edited_gff_to_fasta_with_name_blast_xml','r') as xmlblast2:
     soaptair_record = NCBIXML.parse(xmlblast2)
     species = find_the_subject_db('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAP_ragatag_with_read_correction_to_brassicas_blast/TAIR10_cdna_20101214_updated_SOAP_ragtag_BO_with_read_correction.fasta_xml')
     for record in soaptair_record:
         for description in record.descriptions:
             for alignment in record.alignments:
                 for hsp in alignment.hsps:
                     pre_title  = alignment.title.split(' ')[1]
                     wanted_title = pre_title.split('.')[0]
                     if species == 'arabidopsis':
                         if id_tables['Bnapus'].str.contains(wanted_title).sum():
                             print(wanted_title)
