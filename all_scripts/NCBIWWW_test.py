import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2
import random
'''Rough script to blast the reads that did not align to any of the fastq_screen
databases. Takes a fasta as input (preferably a random sample of a larger fasta)
and needs an input. There are a few hacky bits - mostly the dedup function because
of the for loop in the NCBIXML.parse'''


'''Setting the arguments'''
parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--query', required = True, type = str)
parser.add_argument('--out', required = True, type = str)
args =parser.parse_args()
query = args.query
out = args.out

'''Setting up the output names'''
xml_out = out + '_xml'
summary_out = out + '_xml_summary'
dedup_out = out + '_final'

'''Running the blast - outputting into xml file'''
results_out = open(xml_out,'w+')
fa = SeqIO.parse(open(query,'r'),'fasta')
for seq in fa:
    results_handle = NCBIWWW.qblast('blastn','nt',seq.format('fasta'))
    results_out.write(results_handle.read())


'''This can be changed if needed - for now just outputting what the hits are and
nothing else'''
def xml_summary(xml,out):
    with open(xml,'r') as xmlblast:
        record = NCBIXML.parse(xmlblast)
        descriptions = []
        for rec in record:
            for description in rec.descriptions:
                for alignment in rec.alignments:
                    for hsp in alignment.hsps:
                        out.write(alignment.title + '\n')

'''Because of the for loops in the xml parser, there are often lots of duplicated
lines, this gets rid of them'''
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



xml_summary(xml_out,summary_out)
remove_duplicate_lines(xml_out, dedup_out)
