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
import subprocess

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

results = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/brapa_trichome_to_sinapis_mrna_blast','r')
orfs = '/Users/josephbeegan/Desktop/Work/Sinapis_alba/Test_orfs_fixed.fa'
edit = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/edit_distance_betweeen_sinapis_codingseq_brap_trichomes.txt','r')
trichome_genes = '/Users/josephbeegan/Desktop/Work/Sinapis_alba/brassica_gene_sequences.fa'
sinapis_sequence = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/contigs.fasta','r')
zhang_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/mrna_blast_hits_to_spadesrnaassembly','r')
trichome_names_tair = ['AT5G53200', 'AT3G27920','AT1G79840','AT2G46410']
gene_IDs = []
distance = []
for line in edit:
    gene_IDs.append((line.split()[0] + '_' + line.split()[2]))
    distance.append((int(line.split()[7])))

dictionary = {'Gene_IDs':gene_IDs,'Distance':distance}
blast_df = pd.DataFrame(dictionary)
print(blast_df.head())
order_blast_df = blast_df.sort_values(by=['Distance'])
print(order_blast_df)

#order_blast_df.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/simplified_ordered_sinapis_distance_to_trichomes.txt', sep = '\t', index = False)
trichome_IDs = []
trichome_seq = []
for seq in SeqIO.parse(trichome_genes, 'fasta'):
    trichome_IDs.append(str(seq.id))
    trichome_seq.append(str(seq.seq))

trichome_dict = {'Gene_IDs':trichome_IDs,'Sequence':trichome_seq}
trichome_df = pd.DataFrame(trichome_dict)
print(trichome_df.head())


#Finding the mrna blast hits that match trichome tair IDs in the rnaSPAdes assembly.
rna_assembly_sequences_trichomes = []
best_identity_match = 0
zhang_record = NCBIXML.parse(zhang_blast)
for blast_record in zhang_record:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(alignment.title)
            print(hsp.sbjct)
            rna_assembly_sequences_trichomes.append(hsp.sbjct)
            print(hsp.match)
            print(hsp.query)
            print(hsp.identities)
            print(hsp.positives)
            print('score= ' + str(hsp.score))
            if hsp.positives > best_identity_match:
                best_identity_match = hsp.positives
print('best identity match= ' + str(best_identity_match))



for i in rna_assembly_sequences_trichomes:
    print(i)
#Need to split the fasta file of the rnaSPAdes assembly into seq.seq format so
#That I can find the subject matches in it
rna_no_newline = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdesRNA_assembly_no_newline.fasta','w+')
for seq in SeqIO.parse('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdesRNA_assembly.fasta','fasta'):
    rna_no_newline.write(str(seq.id) + '\n')
    rna_no_newline.write(str(seq.seq) +'\n')




brapa_ID = []
sinapis_seq = []
length = []
records=NCBIXML.parse(results)
outfile = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/brapa_to_sinapis_mrna_summary', 'w+')
for blast_record in records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            for i in trichome_df['Sequence']:
                if levenshtein_ratio_and_distance(hsp.query, i,ratio_calc = True) > 0.70:
                    outfile.write(alignment.title + '\n')
                    outfile.write('Query, match, subject' + '\n')
                    outfile.write(str(trichome_df.loc[trichome_df['Sequence'] == i, 'Gene_IDs'].values) + '\n')
                    brapa_ID.append(trichome_df[trichome_df['Sequence']==i]['Gene_IDs'].values)
                    outfile.write('Length: ' + str(alignment.length) + '\n')
                    outfile.write(hsp.query + '\n')
                    outfile.write(hsp.match + '\n')
                    outfile.write(hsp.sbjct + '\n')
                    sinapis_seq.append(hsp.sbjct)
                    outfile.write('E-value: ' + str(hsp.expect) + '\n')
                    length.append(alignment.length)
outfile.close()

sinapis_dict = {'Gene_IDs':brapa_ID,'Sequence':sinapis_seq,'Length':length}
sinapis_df = pd.DataFrame(sinapis_dict)
print(sinapis_df.head())
print(sinapis_df)
seen = []
for string in sinapis_df['Sequence']:
    string2 = (string[0:15])
    length = list(sinapis_df[sinapis_df['Sequence']==string]['Length'].values)
    print(length[0])
    if string2 not in seen:
        seen.append(string2)
        print(string2)
    else:
        continue
for seq in SeqIO.parse(sinapis_sequence,'fasta'):
    for string in sinapis_df['Sequence']:
        string2 = string[0:15]
        length_wanted = list(sinapis_df[sinapis_df['Sequence']==string]['Length'].values)
        if string2 in seq.seq:
            print(str(sinapis_df[sinapis_df['Sequence']==string]['Gene_IDs']))
            print(str(seq.seq.split(string)[0][0:int(length_wanted[0])]))
#sinapis_df.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/sinapis_sequences_with_best_trichome_matches', sep = '\t', index = False

#print(len(hits))
#missed = set()
#qdict = SeqIO.index(orfs,'fasta')
#print(len(qdict))
#missed = set(qdict.keys()) - set(hits)
#print(len(missed))
#orphan_records = [qdict[name] for name in missed]
#print(orphan_records)
#print(len(qdict))
#misses_prep = qdict[qdict not in hits]
#print(misses_prep)
#misses.add(misses_prep)


#print(misses)
#orphan_records = [q_dict[name] for name in misses]
#print(orphan_records)
for query in SearchIO(blast_list,'blast-text'):
    print(query.query)
#for seq in SeqIO.parse(orfs,'fasta'):
#    print(seq)


for record in NCBIXML.parse(blast_list):
    print(record)





















fruits = {'banana':3,'apple':2, 'mango':1, 'kiwi':5}
cDNA = "/Users/josephbeegan/Downloads/Arabidopsis_thaliana.TAIR10.cdna.all.fa"
CDS = "/Users/josephbeegan/Downloads/Arabidopsis_thaliana.TAIR10.cds.all.fa"
C2H2_list = "/Users/josephbeegan/Desktop/Work/ChIPSeqKNU/C2H2_family_TAIR_IDs_ten.csv"
aminoacids = ['I','M','T','N','F','L','G','S','Y','R','P','W','H','Q','R','A','V','D','E','K']
codon_table = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
conservative_amino_acids = {'Aliphatic':["G","A","V","L","I"],"Cyclic":"P", "Basic":["H","K","R"], "Aromatic":["F","Y","W"],"Acidic":["D","E","N","Q"],"Hydroxyl":["S","C","U","M","T"]}
conservative_amino_acids.get("A")


fix_wanted = []
with open(C2H2_list) as id_name:
    #Just as an aside, the .split() method split a string into a list. Inside the
    #bracekts is the (seperator, maximum number of splits).
    #What this means for below is that it will only take a column of input IDs,
    #any strings outside of that column will be ignored
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_name)
    print(len(wanted))
    for line in wanted:
        fix_wanted.append(line.upper())
#def extract_sequence_by_ID(input,ID_list):
CDSparse = SeqIO.parse(CDS, "fasta")
print(fix_wanted)
addon = '.1'
fix_wanted2 = set([x + addon for x in fix_wanted])
list_wanted = list(fix_wanted2)
print((fix_wanted2))
print(list_wanted[1])
CDS_C2H2 = (seq for seq in SeqIO.parse(CDS, "fasta") if seq.id in fix_wanted2[1])
CDNA_C2H2 = (seq for seq in SeqIO.parse(cDNA, 'fasta')if seq.id in fix_wanted2[1])
#(r for r in SeqIO.parse(CDS, "fasta") if r.id in fix_wanted2)
#for i in SeqIO.parse(records,'fasta'):
#    print(i.id)
#count = SeqIO.write(records, "/Users/josephbeegan/Desktop/TESTTTTT.txt", "fasta")
#trying to find the UTRs by essentially subtracting the CDS from cDNA
#for test in SeqIO.parse(cDNA, 'fasta'):
#    if test.id == 'AT1G55110.1':
#        print(test.seq)
for sequence in SeqIO.parse(CDS,'fasta'):
    if sequence.id == 'AT3G23130.1':
        for seq in SeqIO.parse(cDNA,'fasta'):
            if seq.id == 'AT3G23130.1':
                #out.write("cDNA     " + str(seq.seq))
                oop = seq.seq
                threeprime = seq.seq.split(sequence.seq[0:10])[0]
                fiveprime = seq.seq.split(sequence.seq[-10:])[-1]
                print('3primeUTR of ' + str(seq.id) + ' ' + str(fiveprime) + "\n")
                print('5primeUTR of ' + str(seq.id) + ' ' + str(threeprime) + "\n" + "\n")






fiverprimeUTR = []
IDS = []
outfile = "/Users/josephbeegan/Desktop/UTRfinder.txt"
with open(outfile, 'w') as out:
    for sequence in SeqIO.parse(CDS,'fasta'):
        if sequence.id in list_wanted:
            for seq in SeqIO.parse(cDNA,'fasta'):
                if seq.id in list_wanted:
                    if seq.id == sequence.id:
                    #out.write("cDNA     " + str(seq.seq))
                        oop = seq.seq
                        threeprime = seq.seq.split(sequence.seq[0:10])[0]
                        fiveprime = seq.seq.split(sequence.seq[-10:])[-1]
                        fiverprimeUTR.append(str(fiveprime))
                        IDS.append(seq.id)
                        out.write('5primeUTR of ' + str(seq.id) + ' ' + str(fiveprime) + "\n" + "\n")
                        out.write('3primeUTR of ' + str(seq.id) + ' ' + str(threeprime) + "\n" + "\n")
    out.close()







nuc5uorfs= open("/Users/josephbeegan/Desktop/nucleotide_5prime_uORFs.txt","w+")
numATG = 0


fiveprimeuORFS_lines = open("/Users/josephbeegan/Desktop/fiveprimeuORFS_lines.txt","w+")


print(fiverprimeUTR)
fiveprimeuORFS = []
_dict = dict()
for id,orf in zip(IDS,fiverprimeUTR):
    print("Writing " + str(id) + " to file." )
    fiveprimeuORFS_lines.write(">" + id + "\n")
    key = id
    for eep in orf.split("ATG")[1:]:
        if len(eep) > 15:
                print("ATG" + eep + '\n'+ '\n')
                fiveprimeuORFS_lines.write("ATG" + eep + "\n")
                fiveprimeuORFS.append("ATG" + eep)
                _dict.setdefault(key,[])
                _dict[id].append(eep)
fiveprimeuORFS_lines.close()
print(_dict)
