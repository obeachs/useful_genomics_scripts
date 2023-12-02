import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
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
input_file = "/Users/josephbeegan/Downloads/Arabidopsis_thaliana.TAIR10.cdna.all.fa"
outfile = "/Users/josephbeegan/Desktop/biopytonout.txt"
wanted_ids = "/Users/josephbeegan/Desktop/wanted_ids.txt"




rna_assembly_sequences_trichomes = []
zhang_record = NCBIXML.parse(open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/mrna_blast_hits_to_spadesrnaassembly','r'))
for blast_record in zhang_record:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(alignment.title)
            print(hsp.sbjct)
            rna_assembly_sequences_trichomes.append(hsp.sbjct)
            print(hsp.match)
            print(hsp.query)
#Need to split the fasta file of the rnaSPAdes assembly into seq.seq format so
#That I can find the subject matches in it
for seq in SeqIO.parse('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdesRNA_assembly.fasta','fasta'):
    for subject in rna_assembly_sequences_trichomes:
        print(try_all_matches(seq.seq, subject, 900))









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

for i in SeqIO.parse(count, 'fasta'):
    print(i.id)
#print("Saved %i records from %s to %s" % (count, input_file, outfile))
#for r in SeqIO.parse(input_file, "fasta"):
if count < len(wanted):
    print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))
print(len(wanted))
SEQ1 = []
for seq in SeqIO.parse(outfile, "fasta"):
    print(seq.seq)
    SEQ1.append(seq.seq)
print(len((SEQ1[1])))
try_all_matches(SEQ1[0],SEQ1[1], 90)
