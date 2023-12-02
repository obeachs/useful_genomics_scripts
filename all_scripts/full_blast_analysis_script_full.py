import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2

def auto_blast(blastmode,queryfile, database, outname):
    blast_command = [blastmode, '-query', queryfile, '-db',  database, '-out', outname,  '-num_threads', '8', '-outfmt', '5', '-evalue', '10', '-max_hsps', '100']
    subprocess.run(blast_command)



'''Script to go from gene IDs to BLAST results in a given genome and database'''
'''Needed variables are a list of gene IDS, the cDNA file from which those gene IDS come from,
an output directory and the blast database that they are linked to'''
'''If a blast database needs to be created, use the function makeblastdb'''
'''Default mode is blastn for nucleotides, change blast option to p or x for blastp and blastx
respectively '''

'''Setting up the different variables'''
parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--queryIDS', required = False, default = 'All', type = str)
parser.add_argument('--cdnafile',required = True, default = None, type = str)
parser.add_argument('--blastdatabase', required = True, default = None, type= str)
parser.add_argument('--outfolder', required = True, default = None, type = str)
parser.add_argument('--blast', required = False, default = 'n', type = str)
args =parser.parse_args()



'''Naming of the initial files is important - the script will automatically create
the results and analysis files based off their original names. So for the list of geneIDs
you could have something along the lines of 'WT_MOCK_v_BrPRT6_Mock.txt' '''
'''If the queryIDS parameter is empty all of the ids in the cdna_file will be used'''


blast = args.blast
query = args.queryIDS
'''Important that the cDNA file input into the script matches the input IDs'''
cdna = args.cdnafile
db = args.blastdatabase
out = args.outfolder
if str(out)[-1] != '/':
    out = str(out) + '/'
if blast =='n':
    blastmode = 'blastn'
if blast =='p':
    blastmode = 'blastp'
if blast =='x':
    blastmode = 'blastx'

'''Function that takes a text file list of gene ids and converts them to a python list'''
def excel_convert_to_text(excel):
    outfile = open(id_outname,'w+')
    df = pd.read_excel(excel)
    first_column = df.iloc[:, 0]
    id_list = first_column.tolist()
    for id in id_list:
        outfile.write(str(id)+'\n')




def extract_IDS(file):
    id_list = []
    id_file = open(file,'r')
    for line in id_file:
        id_list.append(line[:-1])
    id_file.close()
    return id_list

def extract_IDS_fasta(fasta):
    id_list = []
    with open(fasta,'r'):
        for seq in SeqIO.parse(fasta,'fasta'):
            id_list.append(seq.description)
    return id_list

'''Function that references the cdna file with the list of ids and generate a new
fasta file with those ids and genes alone'''
print('HELLO HELLO')
print(out)
def generate_fasta(id_list, cdnafile, outfile):
    out = open(outfile, 'w+')
    count = 0
    with open(cdna) as genomefasta:
        count = 0
        for seq in SeqIO.parse(genomefasta,'fasta'):
            ID = seq.description
            if str(ID) in id_list:
                count += 1
                out.write('>' + str(seq.description) + '\n' + str(seq.seq) + '\n')
    out.close()
    if count == 0:
        print('No IDs found matching in cDNA file - closing')
        exit()

'''Function that uses 'subprocess' to make use of the command line to use blastn to blast
the new fasta file to the database given'''

def auto_blast(blastmode,queryfile, database, outname):
    blast_command = [blastmode, '-query', queryfile, '-db',  database, '-out', outname,  '-num_threads', '8', '-outfmt', '5', '-evalue', '10', '-max_hsps', '100']
    subprocess.run(blast_command)

def blast_xml_analysis(id_list, xml,outname):
    out = open(outname,'w+')
    out.write('sinapis_id' + '\t' + 'tair_id' + '\t' + 'sinapis_start' + '\t' + 'sinapis_end' + '\t' + 'hit_length' + '\t' + 'tair_start' + '\t'+'tair_end' + '\n')
    with open(xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        tair_id =alignment.title.split(' ')[0]
                        barley_id = record.query
                        end = hsp.query_start + alignment.length
                        out.write(barley_id + '\t' + tair_id + '\t' + str(hsp.query_start) + '\t' + str(end)  + '\t' + str(alignment.length) + str(hsp.sbjct_start)  + '\t' + str(int(hsp.sbjct_start) + int(alignment.length)) + '\n' )
    out.close()

def make_matches_fasta(id_list, xml,outname):
    out = open(outname,'w+')
    out.write('sinapis_id' + '\t' + 'tair_id' + '\t' + 'sinapis_start' + '\t' + 'sinapis_end' + '\t' + 'hit_length' + '\t' + 'tair_start' + '\n')
    with open(xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        tair_id =alignment.title.split(' ')[0]
                        barley_id = record.query
                        end = hsp.query_start + alignment.length
                        out.write(barley_id + '\t' + tair_id + '\t' + str(hsp.query_start) + '\t' + str(end)  + '\t' + str(alignment.length) + str(hsp.sbjct_start)  + '\n' )
    out.close()


'''Because of the nature of the outputfile, we have to loop through it several times,
resulting in an output file with many duplicate lines. This function will very quickly
sort through the input and produce ann output with no duplicates'''

def remove_duplicate_lines(input,output):
    print('Processing: ' + input)
    input_file = open(input,'r')
    lines_seen = set() # holds lines already seen
    with open(output, "w+") as output_file:
        for each_line in input_file:
            if each_line not in lines_seen: # check if line is not duplicate
                output_file.write(each_line)
                lines_seen.add(each_line)
    input_file.close()

'''Making the output file names'''
if query != 'All':
    queryname_root = query.split('/')[-1]
    queryname_root = queryname_root.split('.')[0]
else:
    queryname_root = cdna.split('/')[-1]+ '_all_ids_BLAST'

id_outname = out +queryname_root +'.txt'
xml_outname = out + queryname_root + '_blast_to_' + db.split('/')[-1] + '_xml'
fasta_outname = out + queryname_root + '_wanted_genes.fa'
matches_outname = xml_outname[0:-4] + '_matches_list.txt'
dedup_outname  = matches_outname.split('.')[0] + '_dedup.txt'

print('HELP ME')
print(fasta_outname)
print('PLS')
print(queryname_root)
print('YARP')
print(query)
print('MAYBE')
print(cdna)
print(cdna.split('/')[-1])

'''Now running the analysis in order'''
'''#1 isolating the ids into a list and generating a fasta file'''
print('Generating fasta file using  '  + query + '  and  ' + cdna)
if query != 'All':
    if query.split('.')[-1] == 'xlsx':
        excel_convert_to_text(query)
        ids = extract_IDS(id_outname)
        generate_fasta(ids, cdna, fasta_outname)
    else:
        ids = extract_IDS(query)
        print(ids)
        generate_fasta(ids, cdna, fasta_outname)
else:
    ids = extract_IDS_fasta(cdna)
    generate_fasta(ids,cdna,fasta_outname)

'''#2 BLAST the new fasta file against the given database to generate xml file'''
print('BLASTING  ' + fasta_outname + '  against  ' + db)
auto_blast(blastmode,fasta_outname, db, xml_outname)


'''#3 Extracting the hits and input IDs from the blast xml file'''
print('ANALYSING  ' + xml_outname)
blast_xml_analysis(ids, xml_outname, matches_outname)



'''#4 Deduplicating lines from the matches file'''
remove_duplicate_lines(matches_outname, dedup_outname)
