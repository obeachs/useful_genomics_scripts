import sys
import Bio
import re
from Bio import SeqIO
from Bio.Seq import Seq
import os
from Bio.Seq import UnknownSeq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import itertools
import argparse
import subprocess

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--exp1', required = True, default = None, type = str)
parser.add_argument('--exp2',required = False, default = None, type = str)
parser.add_argument('--inp1',required = False, default = None, type = str)
parser.add_argument('--inp2',required = False, default = None, type = str)
parser.add_argument('--ref', required = True, default = None, type= str)
parser.add_argument('--outfolder', required = True, default = None, type = str)
args =parser.parse_args()


def get_rootname(string):
    location = string.find('rep')
    newstring = string[0:location+4]
    newstring = newstring.split('/')[-1]
    return newstring


def check_zip(readfile):
    if readfile[-2:] == 'gz':
        print(str(readfile) + ' is zipped, please unzip all read files and re-run the script')
        quit()

    
        
        
inp1 = args.inp1
inp2 = args.inp2
exp1 = args.exp1
exp2 = args.exp2

ref = args.ref

full_list = [inp1,inp2,exp1,exp2]
input_list = [inp1,inp2]
exp_list = [exp1,exp2]

'''Quick check to make sure all files are unzipped'''
# for readfile in full_list:
#     check_zip(readfile)

exp_rootname = get_rootname(exp1)
inp_rootname = get_rootname(inp1)


print(exp_rootname)



'''Setting up the output folders'''
parent_dir = args.outfolder
fastqc_path = os.path.join(parent_dir, 'fastqc_results')
if os.path.isdir(fastqc_path) == False:
    os.mkdir(fastqc_path)
alignment_path = os.path.join(parent_dir, 'aligment_results')
if os.path.isdir(alignment_path) == False:
    os.mkdir(alignment_path)
stats_path = os.path.join(parent_dir, 'alignment_stats')
if os.path.isdir(stats_path) == False:
    os.mkdir(stats_path)
bigwig_path = os.path.join(parent_dir, 'bigwig_results')
if os.path.isdir(bigwig_path) == False:
    os.mkdir(bigwig_path)
peaks_path = os.path.join(parent_dir, 'peaks')
if os.path.isdir(peaks_path) == False:
    os.mkdir(peaks_path)

'''Setting up the commands'''
    
'''COMMAND LIST IF NEEEDED

VARIABLES BEING CALLED BEFORE THEY'RE MADE BREAKS SCRIPT, HAVE ALL COMMANDS HERE FOR QUICK REFERENCE

    
fastqc_command = ['fastqc', '-o', fastqc_path, '-f', reads]
hisat_alignments_command = ['hisat2', '-p', '10', input1, '-x', '~/refs/Hisat2/tair10', '-U', '$R4.fq', '-S', '$R4_hisat2.sam']
bowtie_alignments_command = ['bwa', 'mem', '-t', '2', input1, input2, '>' outsam]
bam_command = ['samtools', 'view' ,'-Sb',outsam, '>', outbam]
sort_command = ['samtools','sort',outbam,outsortedbam]
flagstat_command = ['samtools', 'flagstat', outsortedbam, outflagstat]
rmdup_command = ["samtools",'rmdup','-s',outsortedbam, outrmdup]
begraph_command = ['bedtools', 'genomecov', '-ibam', outrmdup, '-g', ref, '-bg', '>', outbegraph]
'''


'''Aligning the reads based on whether or not it is single end or paired end sequencing'''

print(fastqc_path)
print(alignment_path)
print(stats_path)
print(bigwig_path)

print(ref)


'''Carrying out fastqc on all of the reads, both input and exp'''
# for reads in full_list:
#     fastqc_command = ['fastqc', '-o', fastqc_path, '-f',fastq, reads]
#     subprocess.run(fastqc_command)
    
'''Aligning the reads to the reference genome'''
inp_outsam = os.path.join(alignment_path, inp_rootname + '.sam')
exp_outsam = os.path.join(alignment_path, exp_rootname + '.sam')
inp_bowtie_alignments_command = ['bowtie2', '-p', '4', '-x','/Volumes/sesame/joerecovery/genomes/TAIR10_chr_all','-1',inp1,'-2',inp2, '-S',inp_outsam]
exp_bowtie_alignments_command = ['bowtie2', '-p', '4', '-x','/Volumes/sesame/joerecovery/genomes/TAIR10_chr_all','-1',exp1,'-2',exp2, '-S',exp_outsam]
print('Aligning inputs to reference genome:  ' + str(ref))
print('Results will be found in  ' + str(alignment_path))
subprocess.run(inp_bowtie_alignments_command)
print('Aligning tests to reference genome:  ' + str(ref))
print('Results will be found in  ' + str(alignment_path))
subprocess.run(exp_bowtie_alignments_command)


'''Continuing the alignment process'''
inp_outbam = os.path.join(alignment_path, inp_rootname + '.bam')
inp_bam_command = ['samtools', 'view' ,'-Sb',inp_outsam, '>', inp_outbam]
exp_outbam =  os.path.join(alignment_path, exp_rootname + '.bam')
exp_bam_command = ['samtools', 'view' ,'-Sb',exp_outsam, '>', exp_outbam]
subprocess.run(inp_bam_command)
subprocess.run(exp_bam_command)

inp_outsortedbam = os.path.join(alignment_path, inp_rootname + '_sorted.bam')
inp_sort_command = ['samtools','sort',inp_outbam,inp_outsortedbam]
exp_outsortedbam = os.path.join(alignment_path, exp_rootname + '_sorted.bam')
exp_sort_command = ['samtools','sort',exp_outbam,exp_outsortedbam]
subprocess.run(inp_sort_command)
subprocess.run(exp_sort_command)

inp_outrmdup = os.path.join(alignment_path, inp_rootname + '_sorted_rmdup.bam')
inp_rmdup_command = ["samtools",'rmdup','-s',inp_outsortedbam, inp_outrmdup]
exp_outrmdup = os.path.join(alignment_path, exp_rootname + '_sorted_rmdup.bam')
exp_rmdup_command = ["samtools",'rmdup','-s',exp_outsortedbam, exp_outrmdup]
subprocess.run(inp_rmdup_command)
subprocess.run(exp_rmdup_command)

'''Generating the alignment statistics for troubleshooting and bias detection'''
inp_outflagstat = os.path.join(alignment_path, inp_rootname + '_flagstat.txt')
inp_flagstat_command = ['samtools', 'flagstat', inp_outsortedbam, inp_outflagstat]
exp_outflagstat = os.path.join(alignment_path, exp_rootname + '_flagstat.txt')
exp_flagstat_command = ['samtools', 'flagstat', exp_outsortedbam, exp_outflagstat]
subprocess.run(inp_flagstat_command)
subprocess.run(exp_flagstat_command)

'''Peak calling using macs2, if a different peak caller is needed, edit the command'''
macs2_command = ['macs2','callpeak','-t', exp_outsortedbam, '-c', inp_outsortedbam, '-f', 'BAM', '-g', '1.3e+8', '--outdir', peaks_path, '-n', exp_rootname + '_peaks', '2>', peaks_path + exp_rootname+'_macs2.log']
subprocess.run(macs2_command)


