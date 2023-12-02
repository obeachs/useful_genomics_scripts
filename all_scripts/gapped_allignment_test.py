import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
import numpy as np


def fasta_sequences_recover(fasta):
    sequence_list = []
    for seq in SeqIO.parse(fasta,'fasta'):
        seqs = str(seq.seq)
        sequence_list.append(seqs)
    return(sequence_list)


real_query = "ATGTCAGGCTTAGCATGTAATTCATGTAACAAGGAGTTCGAAGACGACGC"
gap_query_test = "ATGGGGGGTCAGGCTTAGCATGTAATTCATGTAACAAGGAGTTCGAAGACGACGC"
test_query = "ATGTCAGGCTTAGCATGTAATTCATGTAACAAGGGGGGGGGGTTCGAAGACGACGC"
brapa_sequence = "/Users/josephbeegan/Desktop/Work/Sinapis_alba/Brapa_Chr3_sample.fasta"
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
                cost = 0
                 # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
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
def score_match(subject,query):
    new_seq = []
    # for each base in the match
    for i in range(0,len(subject)):
        for i in range(0,len(query)):
            subject_base = subject[i]
            query_base = query[i]
            # then adjust the score up or down depending on
            # whether or not they are the same
            if subject_base == query_base:
                new_seq.append(query_base)
            else:
                n = i + 1
                while n < i + 5:
                    new_query_base = query[n]
                    if subject_base == new_query_base:
                            new_seq.append(new_query_base)
                    else:
                        continue
    print(new_seq)
#score_match(real_query, gap_query_test)
def new_alignment_test(subject, query):
    wanted = 0
    final_length = None
    final_subject_start = None
    final_query_start = None
    for j in range(0,len(subject)):
        for i in range(0,len(query)):
            if query[0:i] in subject[0:j]:
                print("trying to find " + query[0:i] + " in subject")
            else:
                n = i +1
                if query[i:n] in subject[0:j]:
                    print("trying to find " + query[i:n] + " in subject")
def print_one_gap_two_gap_three(list):
    string = ''
    for name in list:
        if str(name).isdigit():
            string = string + '_'*int(name)
        else:
            string = string + name
    print(string)


def matrix_alignment_test(subject,query):
    rows = len(subject)+1
    cols = len(query)+1
    distance = np.zeros((rows,cols),dtype = int)
    final_alignment_list = []
    longest_alignment = ''
    second_alignment = ''
    full_alignment = []
    gap = 0
    gap_index = []
    query_end = 0
    subject_end = 0
    # Populate matrix of zeros with the indeces of each character of both strings
    for i in range(1, rows):
        for k in range(1,cols):
            distance[i][0] = i
            distance[0][k] = k
    for col in range(1, cols):
        for row in range(1, rows):
            if subject[0:row-1] == query[0:col-1]:
                aligned = query[0:col-1]

                if len(aligned) > len(longest_alignment):
                    longest_alignment = aligned
                    query_end = col -1
                    subject_end = row - 1
                    full_alignment.append(aligned)
                    gap_index.append(gap)
            else:
                for n in range(len(query)):
                    if subject[subject_end:row -1] == query[query_end + n:col-1]:
                        new_aligned = query[query_end + n:col-1]
                        if len(new_aligned) > len(second_alignment):
                            second_alignment = new_aligned
                            full_alignment.append(new_aligned)
                            gap = (n - query_end +2)
                            gap_index.append(gap)
    aligned_index_gap = dict(zip(full_alignment,gap_index))
    final_aligned = []
    for i in range(len(full_alignment)-1):
        if full_alignment[i] in full_alignment[i+1]:
            pass
        else:
            print(full_alignment[i])
            final_aligned.append(full_alignment[i])
            final_aligned.append(gap_index[i+1])
            final_aligned.append(full_alignment[-1])
    print(aligned_index_gap)
    print(final_aligned)
    print_one_gap_two_gap_three(final_aligned)



matrix_alignment_test(sinapis_alba_subject,brapa_query)
for seq in SeqIO.parse(brapa_sequence, "fasta"):
    matrix_alignment_test(seq.seq,gap_query_test)

test_list = ["A","2","A","3",'A']




def try_all_matches(subject, query, num):
    wanted = 0
    final_length = None
    final_subject_start = None
    final_query_start = None
    for subject_start in range(0,len(subject)):
        for query_start in range(0,len(query)):
            for length in range(0,len(query)):
                if (subject_start + length < len(subject) and query_start + length < len(query)):
                    score = score_match(subject, query, subject_start, query_start, length)
                    # only print a line of output if the score is better than some limit
                    if score >= wanted:
                        ratio =levenshtein_ratio_and_distance(subject[subject_start:len(query)],query[query_start:len(query)], ratio_calc = True)
                        wanted = score
                        if ratio > num :
                            final_length = length
                            final_subject_start = subject_start
                            final_query_start = query_start

                                #best_start = subject_start
                                #best_length = length
                        else:
                            continue
    print("Score is " + str(wanted))


for seq in SeqIO.parse(brapa_sequence, "fasta"):
    try_all_matches(seq.seq,test_query, 0.5)
    print(levenshtein_ratio_and_distance(seq.seq[0:len(real_query)], real_query, ratio_calc = True))
for line in best:
    print(line)
