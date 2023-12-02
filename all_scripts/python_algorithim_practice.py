import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import argparse
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import subprocess
snapfasta = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_snap_prediction.txt','r')
trinity = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/Trinity.fasta.masked','r')
spades = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdesRNA_assembly.fasta','r')
#contigs = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/contigs.fasta','r')
blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/sinapis_augustus_mrna_tairxml','r')
trichome_names_tair = ['AT5G53200', 'AT3G27920','AT1G79840','AT2G46410']
outfile = open('/Users/josephbeegan/Desktop/sinapis_contigs_no_newline.fasta','w+')
soap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_augustus_mrna_blast_TAIR10cDNA.xml','r')
soap_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_augustus.mrna')
spades_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/augustus_SPAdes_genomic_assembly/augustu_contigs_Ns.mrna','r')
spades_snap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_snap_blast_xml','r')
spades_snap_genes = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdes_snap.fasta','r')
soap_snap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_snap_blast_xml','r')
#zhang_spades = open('')
ragtag_spades_genome = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/augustus_predictions/SPAdes_genomic_sinapis_ragtag_scaffolds_augustus3.mrna','r')
ragtag_spades_zhang_genome = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/augustus_predictions/SPAdes_zhang_ragtag_sinapis_scaffolds_augustus3.mrna','r')
ragtag_soap_genome = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/augustus_predictions/SOAP_ragtag_sinapis_scaffolds_augustus3.mrna','r')

ragtag_spades_blast = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_spades_augustus_mrna_tairblast_xml_with_IDs','r')
ragtag_spades_blast_record = NCBIXML.parse(ragtag_spades_blast)
ragtag_zhang_blast = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_spades_zhang_augustus_mrna_tairblast_xml_with_IDs','r')
ragtag_zhang_blast_record = NCBIXML.parse(ragtag_zhang_blast)
ragtag_soap_blast = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_spades_zhang_augustus_mrna_tairblast_xml_with_IDs','r')
ragtag_soap_blast_record = NCBIXML.parse(ragtag_soap_blast)
ragtag_spades_prelim = pd.read_csv('/Users/josephbeegan/Desktop/ragtag_spades_df_prelim.txt', sep = '\t')
ragtag_soap_prelim = pd.read_csv('/Users/josephbeegan/Desktop/ragtag_soap_df_prelim.txt', sep = '\t')
ragtag_zhang_prelim = pd.read_csv('/Users/josephbeegan/Desktop/ragtag_zhang_df_prelim.txt', sep = '\t')
print(ragtag_soap_prelim)
print(ragtag_spades_prelim)


def make_query_dict(dataframe,genome):
    test_dict = dict()
    for seq in SeqIO.parse(genome,'fasta'):
        for i in range(len(dataframe)):
            for que in dataframe['queries'][i]:
                if len(que) > 14:
                    if que in seq.seq.upper():
                        test_dict.setdefault(seq.id,[])
                        test_dict[seq.id].append(que)
    return(test_dict)




record = NCBIXML.parse(blast)
soaprecord = NCBIXML.parse(soap_blast)
soap_snap_record = NCBIXML.parse(soap_snap_blast)
spades_snap_record = NCBIXML.parse(spades_snap_blast)
spades_list = []
'''
for blast_record in record:
    for alignment in blast_record.alignments:
        spades_list.append((alignment.title.split(' ')[1]))

soap_list = []
for soap_record in soaprecord:
    for alignment in soap_record.alignments:
        soap_list.append(alignment.title.split(' ')[1])
count = 0
popop = []
for soap in soap_list:
    for spades in spades_list:
        if soap == spades:
            popop.append(soap)
            count += 1
print(count)
seen = []
for i in popop:
    if i not in seen:
        seen.append(i)
    else:
        continue
print(len(seen))'''
def blast_to_dataframe(BLAST_record):
    titles = []
    starts = []
    lengths = []
    matches = []
    queries = []
    for blast_record in BLAST_record:
        for description in blast_record.descriptions:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    '''From here it is only taking the trichome gene_IDs'''
                    for name in trichome_names_tair:
                        if name in alignment.title:
                            titles.append(alignment.title.split(' ')[1])
                            starts.append(hsp.query_start)

                            lengths.append(alignment.length)
                            queries.append(hsp.query.split('-'))



    pre_dict = {'titles':titles,'lengths':lengths,'starts':starts,'queries':queries}
    datafram = pd.DataFrame(pre_dict)
    return(datafram)


soap_snap_df = blast_to_dataframe(soap_snap_record)
soap_augustus_df = blast_to_dataframe(soaprecord)
spades_snap_df = blast_to_dataframe(spades_snap_record)
spades_augutus_df = blast_to_dataframe(record)
ragtag_spades_df = blast_to_dataframe(ragtag_spades_blast_record)
ragtag_spades_df.to_csv('/Users/josephbeegan/Desktop/ragtag_spades_df_prelim.txt', sep = '\t', index=False)
ragtag_zhang_df = blast_to_dataframe(ragtag_zhang_blast_record)
ragtag_zhang_df.to_csv('/Users/josephbeegan/Desktop/ragtag_zhang_df_prelim.txt', sep = '\t', index=False)
ragtag_soap_df = blast_to_dataframe(ragtag_soap_blast_record)
ragtag_soap_df.to_csv('/Users/josephbeegan/Desktop/ragtag_soap_df_prelim.txt', sep = '\t', index=False)
print(ragtag_soap_df.head())
print(ragtag_soap_df.head())
print(ragtag_soap_df.head())

def make_query_dict(dataframe,genome):
    test_dict = dict()
    for seq in SeqIO.parse(genome,'fasta'):
        for i in range(len(dataframe)):
            for que in dataframe['queries'][i]:
                if len(que) > 14:
                    if que in seq.seq.upper():
                        test_dict.setdefault(seq.id,[])
                        test_dict[seq.id].append(que)
    return(test_dict)


soap_dict = make_query_dict(ragtag_soap_df,ragtag_soap_genome)
spades_dict = make_query_dict(ragtag_spades_df,ragtag_spades_genome)
zhang_dict = make_query_dict(ragtag_zhang_df,ragtag_spades_zhang_genome)

print(soap_dict)
print(spades_dict)
print(zhang_dict)

for seq in SeqIO.parse(spades_snap_genes,'fasta'):
    for i in range(len(datafram)):
        for que in datafram['queries'][i]:
            if len(que) > 14:
                if que in seq.seq.upper():
                    test_dict.setdefault(seq.id,[])
                    test_dict[seq.id].append(que)
print(test_dict)
#print(pd.DataFrame.from_dict(test_dict).head())
with open('/Users/josephbeegan/Desktop/spades_snap_queries_file.txt', 'w+') as file:
    file.write(str(test_dict))
datafram.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/zhang_SPAdes__genomic_assembly_trichome_matches', sep = '\t',index=False)
'''Omitting this block but it does work, just very time-consuming'''




def levenshtein_ratio_and_distance(s, t, ratio_calc = False):
    """ levenshtein_ratio_and_distance:
        Calculates levenshtein distance between two strings.
        If ratio_calc = True, the function computes the
        levenshtein distance ratio of similarity between two strings
        For all i and j, distance[i,j] will contain the Levenshtein
        distance between the first i characters of s and the
        first j characters of t
    """
    # Initialize matrix of zeros
    rows = len(s)+1
    cols = len(t)+1
    distance = np.zeros((rows,cols),dtype = int)

    # Populate matrix of zeros with the indeces of each character of both strings
    for i in range(1, rows):
        for k in range(1,cols):
            distance[i][0] = i
            distance[0][k] = k

    # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0 # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
            else:
                # In order to align the results with those of the Python Levenshtein package, if we choose to calculate the ratio
                # the cost of a substitution is 2. If we calculate just distance, then the cost of a substitution is 1.
                if ratio_calc == True:
                    cost = 2
                else:
                    cost = 1
            distance[row][col] = min(distance[row-1][col] + 1,      # Cost of deletions
                                 distance[row][col-1] + 1,          # Cost of insertions
                                 distance[row-1][col-1] + cost)     # Cost of substitutions
    if ratio_calc == True:
        # Computation of the Levenshtein Distance Ratio
        Ratio = ((len(s)+len(t)) - distance[row][col]) / (len(s)+len(t))
        return Ratio
    else:
         print(distance) # Uncomment if you want to see the matrix showing how the algorithm computes the cost of deletions,
         #insertions and/or substitutions
         #This is the minimum number of edits needed to convert string a to string b



#AT1G79840_match = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/AT1G79840_match_region.txt','w+')
#wanted_id = ''
#wanted_score = 0

#score = 0
#for seq in SeqIO.parse(contigs,'fasta'):
#    for i in range(len(seq.seq) - datafram['lengths'][0]):
#        window = seq.seq[i:i+datafram['lengths'][0]]
#        if levenshtein_ratio_and_distance(window,datafram['queries'][0],ratio_calc = True) > score:
#            score = levenshtein_ratio_and_distance(window,datafram['queries'][0],ratio_calc = True)
#            wanted_id = seq.id
#print(wanted_id)
#print(wanted_score)
#AT1G79840_match.write(wanted_id + '\n' + wanted_score)
length = 0
for seq in SeqIO.parse(contigs,'fasta'):
    for j in reversed(range(len(datafram['queries'][0]))):
        if datafram['queries'][0][0:j] in seq.seq:
            if len(datafram['queries'][0][0:j]) > length:
                length = len(datafram['queries'][0][0:j])
                print(datafram['titles'][0])
                print(seq.id)
                print(datafram['queries'][0][0:j])



outfile = open('/Users/josephbeegan/Desktop/sinapis_contigs_no_newline.fasta','w+')
for seq in SeqIO.parse(contigs,'fasta'):
    outfile.write(seq.id)
    print(seq.id)
    outfile.write(seq.seq)
    for match in matches:
        for meme in match:
            if len(meme) > 10:
                if meme in seq.seq:
                    print(seq.id)






print(onehotmatrix_genes('ATGCTGCNNNNN'))
def levenshtein_ratio_and_distance(s, t, ratio_calc = False):
    """ levenshtein_ratio_and_distance:
        Calculates levenshtein distance between two strings.
        If ratio_calc = True, the function computes the
        levenshtein distance ratio of similarity between two strings
        For all i and j, distance[i,j] will contain the Levenshtein
        distance between the first i characters of s and the
        first j characters of t
    """
    # Initialize matrix of zeros
    rows = len(s)+1
    cols = len(t)+1
    distance = np.zeros((rows,cols),dtype = int)

    # Populate matrix of zeros with the indeces of each character of both strings
    for i in range(1, rows):
        for k in range(1,cols):
            distance[i][0] = i
            distance[0][k] = k

    # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0 # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
            else:
                # In order to align the results with those of the Python Levenshtein package, if we choose to calculate the ratio
                # the cost of a substitution is 2. If we calculate just distance, then the cost of a substitution is 1.
                if ratio_calc == True:
                    cost = 2
                else:
                    cost = 1
            distance[row][col] = min(distance[row-1][col] + 1,      # Cost of deletions
                                 distance[row][col-1] + 1,          # Cost of insertions
                                 distance[row-1][col-1] + cost)     # Cost of substitutions
    if ratio_calc == True:
        # Computation of the Levenshtein Distance Ratio
        Ratio = ((len(s)+len(t)) - distance[row][col]) / (len(s)+len(t))
        return Ratio
    else:
        # print(distance) # Uncomment if you want to see the matrix showing how the algorithm computes the cost of deletions,
        # insertions and/or substitutions
        # This is the minimum number of edits needed to convert string a to string b
        return "The strings are {} edits away".format(distance[row][col])


#outfile = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/distance_between_trinity_and_spadesrna.txt','w+')
#for seq in SeqIO.parse(trinity,'fasta'):
#    for seq2 in SeqIO.parse(spades,'fasta'):
#        outfile.write(str(seq.id) + ' vs.' + str(seq2.id) + ' = ' + str(levenshtein_ratio_and_distance(str(seq.seq),str(seq2.seq),ratio_calc=True)))



seq = 'ATGGGATTCCAGTTGACCAAA'
seen = []
count = 0
for i in seq:
    if i not in seen:
        seen.append(i)
        count += 1
code = dict()
for i,j in zip(range(count),seen):
    print(str(j) + ' = ' + str(i+1))
    code.setdefault(j,)
    code[j] = i+1


print(code)
print(code.items())


matrix = np.zeros((len(seq),int(count)),dtype = int)
print(matrix)
for line,base in zip(matrix,seq):
    line[(code.get(base)) - 1 ] = 1
print(matrix)

test_string = 'Add to the PATH on Mac OS X 10.8 Mountain Lion and up'
def onehotmatrix(seq):
    seen = []
    count = 0
    for i in seq:
        if i not in seen:
            seen.append(i)
            count += 1
    code = dict()
    for i,j in zip(range(count),seen):
        print(str(j) + ' = ' + str(i+1))
        code.setdefault(j,)
        code[j] = i+1
    matrix = np.zeros((len(seq),int(count)),dtype = int)

    for line,base in zip(matrix,seq):
        line[(code.get(base)) - 1 ] = 1
    return(matrix)
print(onehotmatrix(test_string))

def onehotmatrixgenes(seq):
    code = {'A':1,'G':2,'T':3,'C':4}
    matrix = np.zeros((len(seq),4)
    print(matrix)
    #for line,base in zip(matrix,seq):
    #    line[(code.get(base)) - 1 ] = 1
    #    if base not in code:
    #        line = line
    return(matrix)

print(onehotmatrix_genes('ATGCTGCNNNNN'))
