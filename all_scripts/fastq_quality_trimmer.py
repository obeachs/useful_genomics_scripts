import sys
import Bio
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
import argparse

pwd = os.getcwd()
parser = argparse.ArgumentParser(description='Remove FASTQ reads below a certain quality.')
parser.add_argument('reads', type=str)
parser.add_argument('q', type = int, default = 30)
parser.add_argument('out', type = str)
parser.add_argument('dir', type=str)
args = parser.parse_args()

num = 0
for record in SeqIO.parse(args.reads,'fastq'):
    num += 1
print("There are %i reads in total " % num)

quality_reads = (
    reads
    for reads in SeqIO.parse(args.reads, 'fastq')
    if min(reads.letter_annotations['phred_quality']) >= args.q
)

out = SeqIO.write(quality_reads,args.dir+args.out+ "_score_"+str(args.q)+".fastq", 'fastq')

print("%i reads passed the check of q score " % out)
