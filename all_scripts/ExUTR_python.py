
'''Trying to remake the perl EXUTR script in python'''
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
import itertools
import argparse
import subprocess
import os

sequence = Seq('AAATTGGGGGGAGGGTTAAAAAAAGCCCCCCCGG')
print(sequence.translate())
## Generate a FASTA file of transcripts with potential ORFs in Swissprot.
## Generate a FASTA file of corresponding ORFs of those transcripts.
## Annotate these transcripts (Swissprot ID will be shown on the sequence title)
'''We have to use record.query to extract the title of the blast query used'''
#blast_file ='/Users/josephbeegan/Desktop/merged_test_blast_xml'
#test_file = '/Users/josephbeegan/Desktop/Work/fake_fasta.txt'
brassica_transcriptome = '/Volumes/seagatedrive/stringtie_brapa_transcripts.fasta'
# print(str(blast_file))
# with open(blast_file) as blast_file2:
#     blast_file_record = NCBIXML.parse(blast_file2)
#     for record in blast_file_record:
#         for description in record.descriptions:
#             for alignment in record.alignments:
#                 for hsp in alignment.hsps:
#                     print(alignment.title.split(' ')[0])
#                     print(record.query)
#                     print(hsp.query)
    #def make_final_transcript_and_orf_files(blast_xml_file, transcript_fasta):
def make_final_files(blast_xml,transcripts_fasta):
    Swissprot_titles = []
    transcript_titles = []
    transcript_orfs = []
    transcript_dna = []
    with open(blast_xml) as blast_file:
            blast_file_record = NCBIXML.parse(blast_file)
            for record in blast_file_record:
                for description in record.descriptions:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            print(alignment.title.split(' ')[0])
                            Swissprot_titles.append(alignment.title.split(' ')[0])
                            print(record.query)
                            transcript_titles.append(record.query)
                            print(hsp.query)
                            transcript_orfs.append(hsp.query)
    print(len(transcript_titles))
    print(transcript_titles)
    with open(transcripts_fasta) as handle:
        for seq in SeqIO.parse(handle,'fasta'):
            if str(seq.id) in transcript_titles:
                transcript_titles.append(str(seq.id))
                transcript_dna.append(str(seq.seq))
    transcript_file_name = str(transcripts_fasta.split('.')[0] + '_with_orfs.fa')
    orf_file_name = str(blast_xml + '_with_orfs')
    orf_output = open(orf_file_name,'w+')
    transcript_output = open(transcript_file_name,'w+')
    print('TRANSCRIPT ORFS')
    print(transcript_orfs)
    print('SWISSPROT TITLES')
    print(Swissprot_titles)
    for x, y in zip (Swissprot_titles,transcript_orfs):
        orf_output.write(">" + x + '\n' + y + '\n')
    for w,z in zip(Swissprot_titles, transcript_dna):
        transcript_output.write(">" + w + '\n' + z + '\n')
    orf_output.close()
    transcript_output.close()
#make_final_files(blast_file,test_file)
#make_final_files(blast_file, test_file, '/Users/josephbeegan/Desktop/function_test_tmp')



blastdb = '/Volumes/seagatedrive/uniprot'
parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--filename', required = False, default = '/Users/josephbeegan/Desktop/', type=str)
parser.add_argument('--threads',required = False, default = 2, type = int)
#parser.add_argument('--threeprimeUTR', required = True, default = None, type=argparse.FileType('r'))
args =parser.parse_args()
# brapa_transcripts = '/Volumes/seagatedrive/stringtie_brapa_transcripts.fasta'
# fart = 'AAATTGGGGGGAGGGTTAAAAAAAGCCCCCCCGG'
# seq = Seq('ATGGGGGAAAA')
# trichome_names_tair = ['AT5G53200', 'AT3G27920','AT1G79840','AT2G46410']
# spades_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/augustus_SPAdes_genomic_assembly/augustu_contigs_Ns.mrna','r')

# for nuq in seq:
#     print(nuq)
# for nuc in [seq,seq.reverse_complement()]:
#     print('nuc is ' + str(nuc))
#     for frame in range(3):
#         trans = str(nuc[frame:].translate())
#         print(trans)
# print(args.filename)
# print(str(test_file))
braseq=Seq('CTTTTGTCGCTATCATCCCCACAAGTCTTGAAAGGGGAGCCTTATGATCAAAAGGAACCCCCAAAAACCCTCTCAAACGCATCGCTTTTGTCTGCCACTTTGAGGATCCCAAACAAACCACACCACACTGGTTGTTCATGGATTCTGGGCTTCAGCATCTCGCATTGCGCTTCTTCATTCTTCTCTGCTGCCTCTTTCAAGCCACCGCCGAGGATGGTAGCTGGAAGATAGCCACAGCCACGCTCTCCAGGGACAAGGACGGCTCCTCCTCTGTCACCACTGGAGGCGCTTGTGGGTATGGGGATCTGAGGCAGAGCAGCTATGGCGGGTACAGCGCAGGCCTGAGCGGGAAGCTGTTCAACAGGGGAAGCAGCTGCGGAGCTTGCCTGGAAGTGCGGTGTGTGAACCACATAAGATGGTGCCTTCAAGGCAGCCCCTCCGTGGTGGTCACAGCCACAGATTTCTGTCCTCCCAATTCAGGCCTTTCCTCCGATTACGGAGGTTGGTGCAACTTCCCAAAGGAGCATTTGGAACTATCTCATGCGGCCTTCACAGGGATTGCAGAAACCAGAGCCGAGATAATACCTGTGCAGTACAGGAGGGTCAAGTGTGGGCGGAGAGGCGGGGTGAGATTCAGCTTGAACGGGAGCTCCCACTTCTTTCAGGTCTTGATAAGCAATGTGGGACTCGATGGGGAAGTGGTGGGAGTGAAAGTGAAGGGGCATACAACGGCTTGGATTCCAATGGCCCGAAACTGGGGACAGAACTGGCACTCCTCACTCGATCTCATCGGACAGTCTCTTTCTTTCGAAGTCACTCTCATCGGCGGCAAAACCATTGCCTCTTACGACGTGGCTCCCCCTTATTGGCGCTTTGGAATGACCTACCAAGGAAAGCAGTTCCTTTCGTGACTACCTCTCCTCTTTTCCTTTCTTTATATAAAGTTATGATTCACTTTCGGTGGTCAAACACAATTTGTGGCCTTTATCATGACTTTTAAGACATACTAGTAATCACTTCAGTTCAAGGTATAAATGACTTTACCATTGTCTTTATATATATAATAATAATAGTAATATCCTAC')
'''.find() method will return the location of the string we want to find in the sequence'''

def extract_queries_without_blast_match(blast_xml,query_fasta):
    q_dict = SeqIO.index(query_fasta,'fasta')
    hits = []
    for record in NCBIXML.parse(open(blast_xml)):
        if record.alignments:
            hits.append(record.query.split()[0])

        misses = set(q_dict.keys()) - set(hits)
        orphan_records = [q_dict[name] for name in misses]
    return misses

def find_orfs_in_transcripts_fasta_self(transcript):
    with open(transcript,'r') as handle:
        for seq in SeqIO.parse(handle, "fasta"):
            for i in range(3):
                translated_seq = str(seq.seq[i:].translate())
                print(translated_seq)
                list_of_potential_orfs = translated_seq.split('*')
    return list_of_potential_orfs
# print(find_orfs_in_transcripts_fasta_self(test_file))


def find_orfs_in_revcom_transcripts_fasta_self(transcript):
    with open(transcript,'r') as handle:
        for seq in SeqIO.parse(handle,'fasta'):
            for i in range(3):
                revcom_translated_seq = str(seq.seq.reverse_complement()[i:].translate())
                print(revcom_translated_seq)
                revcom_list_of_potential_orfs = revcom_translated_seq.split('*')
    return revcom_list_of_potential_orfs
# print(find_orfs_in_revcom_transcripts_fasta_self(test_file))



def find_orfs_in_transcripts_self(sequence):
    for i in range(3):
        translated_seq = str(sequence[i:].translate())
        list_of_potential_orfs = translated_seq.split('*')
    return list_of_potential_orfs

def find_orfs_in_revcom_transcripts_self(sequence):
    for i in range(3):
        translated_seq = str(sequence[i:].translate())
        list_of_potential_orfs = translated_seq.split('*')
    return list_of_potential_orfs
print(find_orfs_in_transcripts_self(Seq('AAATTGGGGGGAGGGTTAAAAAAAGCCCCCCCGG')))
print(find_orfs_in_transcripts_self(braseq))
print(find_orfs_in_revcom_transcripts_self(braseq))
#list = find_orfs_in_transcripts_self(braseq)
#print(list)
#print('Here is the max ' + max(list,key = len) + ' it is ' + str(len(max(list,key = len))) + ' long')

def get_seqids(fasta):
    IDs = []
    with open(fasta,'r') as handle:
        for seq in SeqIO.parse(handle,'fasta'):
            IDs.append(seq.id)
    return IDs


def capture_orfs_ids_in_transcript_file(transcript_fasta):
    seq_ids = get_seqids(transcript_fasta)
    count = 0
    orf_list_forward = []
    orf_list_reverse = []
    final_orf1 = []
    final_orf2 = []
    with open(transcript_fasta,'r') as handle:
        for seq in SeqIO.parse(handle,'fasta'):
            count += 1
            #get_seqids(transcript_fasta)
            orf_list_forward.append(max(find_orfs_in_transcripts_self(seq.seq),key= len))
            orf_list_reverse.append(max(find_orfs_in_revcom_transcripts_self(seq.seq.reverse_complement()),key=len))
    outfile_IDs = open('brapa_IDs','w+')
    seq_orf = open('/Users/Desktop/brapa_orfs', 'w+')
    tmp_orf1 = open('/Users/Desktop/brapa_tmp_orf1.fa','w+')
    tmp_orf2 = open('/Users/Desktop/brapa_tmp_orf2','w+')
    print('forward = ')
    print(orf_list_forward)
    print('reverse = ')
    print(orf_list_reverse)
    for i in range(len(orf_list_forward)):
        if len(orf_list_forward[i]) > len(orf_list_reverse[i]):
            orf1 = orf_list_forward[i]
            orf2 = orf_list_reverse[i]
        else:
            orf2 = orf_list_forward[i]
            orf1 = orf_list_reverse[i]
        final_orf1.append(orf1)
        final_orf2.append(orf2)

    # for forward_line,reverse_line in orf_list_forward,orf_list_reverse:
    #     print('Here is an item in the forward_orf list')
    #     print(forward_line)
    #     print('Here is an item in the reverse_orf list')
    #     print(reverse_line)
    #     longest_forward = max(forward_line,key=len)
    #     longest_revcom = max(reverse_line,key=len)
    #     print(forward_line)
    #     print(reverse_line)
    #     if len(longest_forward) > len(longest_revcom):
    #             orf1 = longest_forward
    #             orf2 = longest_revcom
    #     else:
    #         orf1 = longest_revcom
    #         orf2 = longest_forward
    #     final_orf1.append(orf1)
        final_orf2.append(orf2)
    for i in range(len(final_orf1)):
        tmp_orf1.write(">" + seq_ids[i] + '\n' +  final_orf1[i] +'\n')
        tmp_orf2.write('>' + seq_ids[i] + '\n' + final_orf2[i] + '\n')
    outfile_IDs.close()
    seq_orf.close()
    tmp_orf1.close()
    tmp_orf2.close()

# capture_orfs_ids_in_transcript_file(test_file)


def mergexml(list_of_split_files, output_file):
        """Merging multiple XML files is non-trivial and must be done in subclasses."""
        if len(list_of_split_files) == 1:
            #For one file only, use base class method (move/copy)
            return Text.merge(list_of_split_files, output_file)
        out = open(output_file, "w+")
        h = None
        for f in list_of_split_files:
            h = open(f)
            body = False
            header = h.readline()
            if not header:
                out.close()
                h.close()
                raise ValueError("BLAST XML file %s was empty" % f)
            if header.strip() != '<?xml version="1.0"?>':
                out.write(header) #for diagnosis
                out.close()
                h.close()
                raise ValueError("%s is not an XML file!" % f)
            line = h.readline()
            header += line
            if line.strip() not in ['<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">',
                                    '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">']:
                out.write(header) #for diagnosis
                out.close()
                h.close()
                raise ValueError("%s is not a BLAST XML file!" % f)
            while True:
                line = h.readline()
                if not line:
                    out.write(header) #for diagnosis
                    out.close()
                    h.close()
                    raise ValueError("BLAST XML file %s ended prematurely" % f)
                header += line
                if "<Iteration>" in line:
                    break
                if len(header) > 10000:
                    #Something has gone wrong, don't load too much into memory!
                    #Write what we have to the merged file for diagnostics
                    out.write(header)
                    out.close()
                    h.close()
                    raise ValueError("BLAST XML file %s has too long a header!" % f)
            if "<BlastOutput>" not in header:
                out.close()
                h.close()
                raise ValueError("%s is not a BLAST XML file:\n%s\n..." % (f, header))
            if f == list_of_split_files[0]:
                out.write(header)
                old_header = header
            elif old_header[:300] != header[:300]:
                #Enough to check <BlastOutput_program> and <BlastOutput_version> match
                out.close()
                h.close()
                raise ValueError("BLAST XML headers don't match for %s and %s - have:\n%s\n...\n\nAnd:\n%s\n...\n" \
                                 % (list_of_split_files[0], f, old_header[:300], header[:300]))
            else:
                out.write("    <Iteration>\n")
            for line in h:
                if "</BlastOutput_iterations>" in line:
                    break
                #TODO - Increment <Iteration_iter-num> and if required automatic query names
                #like <Iteration_query-ID>Query_3</Iteration_query-ID> to be increasing?
                out.write(line)
            h.close()
        out.write("  </BlastOutput_iterations>\n")
        out.write("</BlastOutput>\n")
        out.close()


capture_orfs_ids_in_transcript_file(brassica_transcriptome)
'''Probably need to figure out the correct inputs for the blast query, need to figure out what the best way to do this is
If I separate the longest orfs of each direction and put them into a fasta file with the initial transcript fasta'''


blast_command = ['blastp', '-query', '/Users/josephbeegan/Desktop/test.txttmp_orf1.fa', '-db', '/Volumes/seagatedrive/uniprot_sprot.fasta', '-out', '/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf1_xml', '-num_descriptions', '1', '-num_alignments', '1', '-num_threads','1', '-evalue', '0.000001', '-outfmt', '5']
subprocess.run(blast_command)

    #os.system('blastp -query ~/Desktop/test.txttmp_orf1.fa -db /Volumes/seagatedrive/uniprot_sprot.fasta -out /Users/josephbeegan/Desktop/blast_text_orf1  -num_descriptions 1 -num_alignments 1 -evalue 0.0001')
'''Now we need to parse the BLAST output and identify any queries without a BLAST hit'''

list_of_split_files = '/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf1_xml'
unblasted_IDS = (extract_queries_without_blast_match(list_of_split_files,test_file))


'''With the information about which of the longer orfs do not have blast hits, we take the corresponding orf from the shorter orf fasta and blast it against the same database'''
unblasted_orfs = open(args.filename + 'tmp_orf1_unblasted.fa','w+')
with open(args.filename + 'tmp_orf2','r') as handle:
    for seq in SeqIO.parse(handle,'fasta'):
        if seq.id in unblasted_IDS:
            unblasted_orfs.write(">" + str(seq.id) + '\n' + str(seq.seq) + '\n')
    unblasted_orfs.close()
orforf = '/Users/josephbeegan/Desktop/test.txttmp_orf1_unblasted.fa'
new_blast_command = ['blastp', '-query',orforf, '-db', '/Volumes/seagatedrive/uniprot_sprot.fasta', '-out', '/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf2_xml', '-num_descriptions', '1', '-num_alignments', '1', '-num_threads','1', '-evalue', '0.000001', '-outfmt', '5']
subprocess.run(new_blast_command)




blast1 = '/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf1_xml'
blast2 = '/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf2_xml'
blast_list = ['/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf1_xml','/Users/josephbeegan/Desktop/test.txttmp_orf1_tmp_blast_orf2_xml']

mergexml(blast_list, 'merged_brapa_xml')

## Generate a FASTA file of transcripts with potential ORFs in Swissprot.
## Generate a FASTA file of corresponding ORFs of those transcripts.
## Annotate these transcripts (Swissprot ID will be shown on the sequence title)
blast_file = open('/Users/josephbeegan/Desktop/merged_test_blast_xml','r')
blast_file_record = NCBIXML.parse(blast_file)
for record in blast_file_record:
    for description in blast_record.descriptions:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print(hsp.title)








forward_orfs = find_orfs_in_transcripts_self(braseq)
revcom_orfs = find_orfs_in_transcripts_self(braseq.reverse_complement())
all_braseq_orf = forward_orfs + revcom_orfs
print(all_braseq_orf)
print(max(all_braseq_orf))





'''39450 read ids within this transcript file'''
transcripts = open('/Volumes/seagatedrive/stringtie_brapa_transcripts.fasta','r')
transcript_ids = []
count = 0
for sequ in SeqIO.parse(transcripts,'fasta'):
    if str(sequ.id) in '>STRG.1.2':
        find_orfs_in_transcripts(sequ.seq)



#parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
#parser.add_argument('--peakfileone', required = True, default = None, type = str, help = 'Homer-annotated ChIP-seq peak file')
#parser.add_argument('--peakfiletwo', required = True, default = None, type= str, help = 'Homer-annotated ChIP-seq peak file')
#parser.add_argument('--outdir', required = False, default = None, type=str, help = 'Output directory - will default to current working directory if not specified')
#parser.add_argument('--gene_family_info', required = False, default = None, type=str, help ='Type anything after this operator to request info about gene family distribution in the shared peaks')
#args =parser.parse_args()


#reformat_fasta_file = `sed 's/>/>lcl|/' $file_input > tmp_rearrange.fa`
#subprocess.Popen(reformat_fasta_file)

orfall = open('>tmp_all_transcripts_id','w+')
orflong = open('>tmp_orf1','w+')
orfshort = open('>tmp_orf2','w+')

#for sequence in SeqIO.parse(args.inputfasta,'fasta'):
#    orfall.write(sequence.id + '\n')
#    forward = seq.translate()
#    reverse = seq.reverse_complement()
#    reverse = reverse.translate()
#    if len(forward) > len(reverse):
#        orflong.write(">" + seq.id + '\n' + forward)
#        orfshort.write(">" + seq.id + '\n' + reverse)
#    else:
#        orflong.write(">" + seq.id + '\n' + reverse)
#        orfshort.write(">" + seq.id + '\n' + forward)

brapa_hairy_genes = "/Users/josephbeegan/Desktop/Work/Sinapis_alba/brassica_gene_sequences.fa"
