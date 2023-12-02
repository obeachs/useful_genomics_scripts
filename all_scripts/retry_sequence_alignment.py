import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
queryA = 'actgatcgattgatcgatcgatcg'
subjectA = 'tttagatcgatctttgatc'
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
            score = score - 1
    return score

def pretty_print_match(subject, query, subject_start, query_start, length):

    # first print the start/stop positions for the subject sequence
    print(str(subject_start) + (' ' * length) + str(subject_start+length))

    # then print the bit of the subject that matches
    print(' ' + subject[subject_start:subject_start+length])

    # then print the bit of the query that matches
    print(' ' + query[query_start:query_start+length])

    # finally print the start/stop positions for the query
    print(str(query_start) + (' ' * length) + str(query_start+length))

    print('\n--------------------\n')

def try_all_matches(subject, query, score_limit):
    for subject_start in range(0,len(subject)):
        for query_start in range(0,len(query)):
            for length in range(0,len(query)):
                if (subject_start + length < len(subject) and query_start + length < len(query)):
                    score = score_match(subject, query, subject_start, query_start, length)
                    # only print a line of output if the score is better than some limit
                    if (score >= score_limit):
                        print('Score : ' + str(score))
                        pretty_print_match(subject, query, subject_start, query_start, length)
input_file = "/Users/josephbeegan/Downloads/Arabidopsis_thaliana.TAIR10.cdna.all.fa"
outfile = "/Users/josephbeegan/Desktop/biopytonout.txt"
wanted_ids = "/Users/josephbeegan/Desktop/wanted_ids.txt"

with open(wanted_ids) as id_name:
    #Just as an aside, the .split() method split a string into a list. Inside the
    #bracekts is the (seperator, maximum number of splits).
    #What this means for below is that it will only take a column of input IDs,
    #any strings outside of that column will be ignored
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_name)
    print(wanted)
print("Found %i unique identifiers in %s" % (len(wanted), wanted_ids))

records = (r for r in SeqIO.parse(input_file, "fasta") if r.id in wanted)
count = SeqIO.write(records, outfile, "fasta")
#print("Saved %i records from %s to %s" % (count, input_file, outfile))
#for r in SeqIO.parse(input_file, "fasta"):
if count < len(wanted):
    print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))
print(len(wanted))
SEQ1 = []
for seq in SeqIO.parse(outfile, "fasta"):
    print(seq.seq)
    SEQ1.append(seq.seq)
print(SEQ1[0])
try_all_matches(SEQ1[0],SEQ1[1], 60)
