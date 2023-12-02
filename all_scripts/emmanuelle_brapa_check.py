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
import matplotlib.pyplot as plt
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2


id_list = []
id_file = open('/Volumes/sesame/joerecovery/emmanuelle_brapa_ids_2/WT_flg22_v_BrPRT6_flg22.txt','r')
for line in id_file:
    id_list.append(line[:-1])
print(id_list)
id_file.close()
genome = '/Volumes/sesame/joerecovery/emmanuelle_blast_results_2021/genome_assemblies_cds_fasta(1)/ncbi-genomes-2021-11-04/Brapa_v3_genes'

'''Preparing the new fasta to be blasted against TAIR to see '''
# required_genes_fasta = open('/Volumes/sesame/joerecovery/emmanuelle_blast_results_2021/WT_flg22_v_BrPRT6_flg22_wanted.fa','w+')



with open(genome) as genomefasta:
    for seq in SeqIO.parse(genomefasta,'fasta'):
        location = seq.description.find('GeneID')
        ID = seq.description[location+7:location+16]
        print(ID)
        if str(ID) in id_list:
            print('yay')
            required_genes_fasta.write('>' + str(seq.description) + '\n' + str(seq.seq) + '\n')


'''Analysing the blast results to extract the matches to TAIR genome'''

outfile = open('/Volumes/sesame/joerecovery/emmanuelle_blast_results_2021/WT_flg22_v_BrPRT6_flg22_wanted_blast_to_TAIR10_cdna_20101214_updated_matches_list.txt','w+')
xml = '/Volumes/sesame/joerecovery/emmanuelle_blast_results_2021/WT_flg22_v_BrPRT6_flg22_wanted_blast_to_TAIR10_cdna_20101214_updated_xml'
with open(xml,'r') as xmlblast2:
        soaptair_record = NCBIXML.parse(xmlblast2)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        brapa_location = record.query.find('GeneID')
                        brapa_id = record.query[brapa_location+7:brapa_location+16]
                        tair_id = alignment.title[0:9]
                        outfile.write(brapa_id + '\t' + tair_id + '\n')