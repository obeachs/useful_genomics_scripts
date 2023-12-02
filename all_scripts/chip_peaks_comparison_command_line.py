import sys
import Bio
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import itertools
import argparse


# parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
# parser.add_argument('--exp1', required = True, default = None, type = str)
# parser.add_argument('--exp2',required = False, default = None, type = str)
# parser.add_argument('--inp1',required = False, default = None, type = str)
# parser.add_argument('--inp2',required = False, default = None, type = str)
# parser.add_argument('--ref', required = True, default = None, type= str)
# parser.add_argument('--outfolder', required = True, default = None, type = str)
# args =parser.parse_args()

bop = 'ref6_alp2-1_input_rep1_1.fq.gz'
def get_rootname(string):
    location = string.find('rep')
    print(string[0:location])
    print(string[0:location+1])
    
get_rootname(bop)



'''Setting up the output folders'''
parent_dir = args.outfolder
fastqc_path = os.path.join(parent_dir, 'fastqc_results')
alignment_path = os.path.join(parent_dir, 'aligment_results')
stats_path = os.path.join(parent_dir, 'alignment_stats')
bigwig_path = os.path.join(parent_dir, 'bigwig_results')


'''Setting up the commands'''
fastqc_command = ['fastqc', '-o', fastqc_path, '-f']
hisat_alignments_command = ['hisat2', '-p', '10', input1, '-x', '~/refs/Hisat2/tair10', '-U', '$R4.fq', '-S', '$R4_hisat2.sam']
bowtie_alignments_command = ['bwa', 'mem', '-t', '2', input1, input2, '>' outsam]
bam_command = ['samtools', 'view' ,'-Sb',outsam, '>', outbam]
sort_command = ['samtools','sort',outbam,outsortedbam]
flagstat_command = ['samtools', 'flagstat', outsortedbam, outflagstat]
rmdup_command = ["samtools",'rmdup','-s',outsortedbam, outrmdup]
begraph_command = ['bedtools', 'genomecov', '-ibam', outrmdup, '-g', ref, '-bg', '>', outbegraph]



'''Aligning the reads based on whether or not it is single end or paired end sequencing'''


