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


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--xmlfolder', required = True, type=str)
parser.add_argument('--output_directory', required = True, type=str)
args = parser.parse_args()

blast_file_folder = args.xmlfolder

directory = args.output_directory

print(str(blast_file_folder))
print(str(directory))

blast_file_list = []
for f in listdir(blast_file_folder):
    blast_file_list.append(blast_file_folder + f)

def find_the_subject_db(xml):
    list = []
    species = ''
    with open(xml,'r') as xmlblast2:
        soaptair_record = NCBIXML.parse(xmlblast2)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        list.append(alignment.title)
    checklist_for_tair = 'AT'
    checklist_for_brapa = 'Bra'
    checklist_for_bo = 'Bo'
    checklist_for_bnapus = 'Bna'

    if 'AT' in list[0]:
        species = 'arabidopsis'
    if 'Bra' in list[0]:
        species = 'brassica_rapa'
    if 'Bo' in list[0]:
        species = 'brassica_oleracea'
    if 'Bna' in list[0]:
        species = 'brassica_napus'

    return species
# print(find_the_subject_db('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAP_ragatag_with_read_correction_to_brassicas_blast/Brassica_napus.AST_PRJEB5043_v1.cds.all.fa_SOAP_ragtag_BO_with_read_correction.fasta_xml'))
# print(find_the_subject_db('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAP_ragatag_with_read_correction_to_brassicas_blast/TAIR10_cdna_20101214_updated_SOAP_ragtag_BO_with_read_correction.fasta_xml'))
trichome_names_tair = ['AT5G53200', 'AT3G27920','AT1G79840','AT2G46410','AT2G30432','AT2G30432']



''' cdna chromosome:Brapa_1.0:A08:350979:352142:-1 gene:Bra030860 gene_biotype:protein_coding transcript_biotype:protein_coding description:AT1G55460 (E=1e-091) | Kin17 DNA-binding protein-related

Make a shortened list based on the ids in the table, then compare the shortened lists/files against each other to
see if the same contigs share homologues across different species'''


id_tables = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/tair_brapa_bnapus_bo_ids.txt',sep='\t')

print(id_tables)

def bnapus_title_identified(string):
    bnapus_location = string.find('Bna')
    wanted_title = string[bnapus_location:bnapus_location+13]
    return wanted_title
def tair_title_identifier(string):
    wanted_title = string.split(' ')[1]
    wanted_title = wanted_title.split('.')[0]
    return wanted_title
def brapa_title_identifier(string):
    brapa_location = string.find('gene:Bra')
    wanted_title = string[brapa_location+5:brapa_location+14]
    return wanted_title
def bo_title_identifier(string):
    bo_location = string.find('gene:Bo')
    wanted_title = string[bo_location+5:bo_location+15]
    return wanted_title


def find_trichome_gene_matches(xmlblast,table):
    filename_shortened = str(xmlblast.split('/')[-1])
    short_summary_file = directory + filename_shortened + '_short_summary.txt'
    long_summary_file = directory + filename_shortened + '_long_summary.txt'

    species = find_the_subject_db(xmlblast)
    print(species)
    print('analysing ' + str(xmlblast))
    print('The species is '+ str(species))
    with open(xmlblast, 'r') as xmlblast2, open(short_summary_file, 'w+') as short_summary, open(long_summary_file,'w+') as long_summary:
        soaptair_record = NCBIXML.parse(xmlblast2)
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if species == 'arabidopsis':
                            wanted_title = tair_title_identifier(alignment.title)
                            if id_tables['TAIR'].str.contains(wanted_title).sum():
                                print(wanted_title)
                                ref_id_list.append(wanted_title)
                                match_id_list.append(record.query)
                                long_summary.write(record.query + ':' +' ' + hsp.query +'\n' + wanted_title + ': ' + hsp.sbjct + '\n' + 'e.value: ' + str(description.e) + '\n' + '\n')
                                short_summary.write(record.query + ':' +' ' + wanted_title +'\n')
                        # if species == 'brassica_rapa':
                        #     wanted_title = brapa_title_identifier(alignment.title)
                        #     if id_tables['Brapa'].str.contains(wanted_title).sum():
                        #         print(wanted_title)
                        #         ref_id_list.append(wanted_title)
                        #         match_id_list.append(record.query)
                        #         # print(wanted_title)
                        #         long_summary.write(record.query + ':' +' ' + hsp.query +'\n' + wanted_title + ': ' + hsp.sbjct + '\n' + 'e.value: ' + str(description.e) + '\n' + '\n')
                        #         # print(wanted_title)
                        #         short_summary.write(record.query + ':' +' ' + wanted_title +'\n')
                        # if species == 'brassica_napus':
                        #     wanted_title = bnapus_title_identified(alignment.title)
                        #     # print('Checking ' + str(xmlblast) + ' against Bnapus IDs')
                        #     if id_tables['Bnapus'].str.contains(wanted_title).sum():
                        #         print(wanted_title)
                        #         ref_id_list.append(wanted_title)
                        #         match_id_list.append(record.query)
                        #         long_summary.write(record.query + ':' +' ' + hsp.query +'\n' + wanted_title + ': ' + hsp.sbjct + '\n' + 'e.value: ' + str(description.e) + '\n' + '\n')
                        #         # print(wanted_title)
                        #         short_summary.write(record.query + ':' +' ' + wanted_title +'\n')
                        # if species == 'brassica_oleracea':
                        #     wanted_title = bo_title_identifier(alignment.title)
                        #     # print('Checking ' + str(xmlblast) + ' against BO IDs')
                        #     if id_tables['BO'].str.contains(wanted_title).sum():
                        #         print(wanted_title)
                        #         ref_id_list.append(wanted_title)
                        #         match_id_list.append(record.query)
                        #         # print(wanted_title)
                        #         long_summary.write(record.query + ':' +' ' + hsp.query +'\n' + wanted_title + ': ' + hsp.sbjct + '\n' + 'e.value: ' + str(description.e) + '\n' + '\n')
                        #         short_summary.write(record.query + ':' +' ' + wanted_title +'\n')

for line in blast_file_list:
    find_trichome_gene_matches(line, id_tables)


summary_file_list = []
for f in listdir(directory):
    if 'short_summary' in f:
        summary_file_list.append(f)
