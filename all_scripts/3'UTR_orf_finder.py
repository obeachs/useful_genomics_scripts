`import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
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
orf_dict = dict()
for id,orf in zip(IDS,fiverprimeUTR):
    print("Writing " + str(id) + " to file." )
    fiveprimeuORFS_lines.write(">" + id + "\n")
    key = id
    for eep in orf.split("ATG")[1:]:
        if len(eep) > 15:
                print("ATG" + eep + '\n'+ '\n')
                fiveprimeuORFS_lines.write("ATG" + eep + "\n")
                fiveprimeuORFS.append("ATG" + eep)
                orf_dict.setdefault(key,[])
                orf_dict[id].append(eep)
fiveprimeuORFS_lines.close()
print(orf_dict)



fiveprimeuORFS_lines = open("/Users/josephbeegan/Desktop/fiveprimeuORFS_lines.txt","r")
fiveprimeuORFS_proteins = []

for uorf in fiveprimeuORFS_lines:
    if ">" in uorf:
        fiveprimeuORFS_proteins.append(uorf)
    else:
        protein= ""
        for i in range(0, len(uorf)-2, 3):
            codon = uorf[i:i+3]
            amino_acid = codon_table.get(codon, '')
            protein = protein + str(amino_acid)
        #print(uorf)
        #print(protein)
        fiveprimeuORFS_proteins.append(protein)

print(fiveprimeuORFS_proteins)




def find_viable_protein_sequences(uORF_proteins, out):
    count = 0
    for prot in uORF_proteins:
        if ">" in prot:
            out.append(prot)
        if prot[1] in aminoacids:
            if prot[2] in aminoacids:
                if prot[3] in aminoacids:
                    viable_uORFS.append(prot)
                    count = count +1
                    out.append(prot)
    print(str(count) + " viable uORFs found")
viable_uORFS = []
find_viable_protein_sequences(fiveprimeuORFS_proteins, viable_uORFS)
print(viable_uORFS)



def score_match(subject, query, subject_start, query_start, length):
    score = 0
    # for each base in the match
    for i in range(0,length):
        # first figure out the matching base from both sequences

        subject_base = subject[subject_start + i]
        query_base = query[query_start + i]
        # then adjust the score up or down depending on
        # whether or not they are the same
        if subject_base == query_base:
            score = score + 1
        else:
            for k in conservative_amino_acids.keys():
                if subject_base in conservative_amino_acids[k]:
                    if query_base in conservative_amino_acids[k]:
                        score = score + 0.5
                    else:
                        score = score -0.5
    return score
def try_all_matches(subject, query, score_limit):
    for subject_start in range(0,len(subject)):
        for query_start in range(0,len(query)):
            for length in range(0,len(query)):
                if (subject_start + length < len(subject) and query_start + length < len(query)):
                    score = score_match(subject, query, subject_start, query_start, length)
                    # only print a line of output if the score is better than some limit
                    if (score >= score_limit):
                        print('Score : ' + str(score))
                        clean_output(subject, query, subject_start, query_start, length)

def clean_output(subject, query, subject_start, query_start, length):

    # first print the start/stop positions for the subject sequence
    print(str(subject_start) + (' ' * length) + str(subject_start+length))

    # then print the bit of the subject that matches
    print(' ' + subject[subject_start:subject_start+length])

    # then print the bit of the query that matches
    print(' ' + query[query_start:query_start+length])

    # finally print the start/stop positions for the query
    print(str(query_start) + (' ' * length) + str(query_start+length))

    print('\n--------------------\n')

reversed_viable_uORFS = viable_uORFS[::-1]

for i,y in zip(viable_uORFS, reversed_viable_uORFS):
    if ">" not in i:
        if ">" not in y:
            if i != y:
                print("Comparing " + i + " to " + y)
                try_all_matches(i,y,4)
