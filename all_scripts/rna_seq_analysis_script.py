import sys
import Bio
import re
import os
from Bio.Seq import UnknownSeq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import itertools
import argparse
import subprocess
import multiprocessing

cores = multiprocessing.cpu_count() # or os.cpu_count()



parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--exp1', required = True, default = None, type = str)
parser.add_argument('--exp2',required = False, default = None, type = str)
parser.add_argument('--ref', required = True, default = None, type= str)
parser.add_argument('--gtf', required = True, default = None, type= str)
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

    
        
        
cont1 = args.cont1
cont2 = args.cont2

exp1 = args.exp1
exp2 = args.exp2

ref = args.ref

stringtie_index = args.gtf


cont_list = [cont1,cont2]
exp_list = [exp1,exp2]
full_list = [cont1,cont2,exp1,exp2]


'''Quick check to make sure all files are unzipped'''
for readfile in full_list:
    check_zip(readfile)

cont_rootname = get_rootname(cont1)
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
stringtie_path = os.path.join(parent_dir, 'stringtie_results')
if os.path.isdir(stringite_path) == False:
    os.mkdir(stringtie_path)

'''Setting up the commands'''
    
'''COMMAND LIST IF NEEEDED

VARIABLES BEING CALLED BEFORE THEY'RE MADE BREAKS SCRIPT, HAVE ALL COMMANDS HERE FOR QUICK REFERENCE




Mostly following the nature protocols method (2016).

If these commands are needed on the command line without the script, just remove all commas and quotation
marks.
 

    
fastqc_command = ['fastqc', '-o', fastqc_path, '-f', reads]
hisat_alignments_command = [' hisat2', '-p', '4', '--dta', '-x',ref, '-1', cont1, '-2', cont2, '-S', outsam]
bowtie_alignments_command = ['bwa', 'mem', '-t', '2', contut1, contut2, '>' outsam]
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


'''Carrying out fastqc on all of the reads, both contut and exp'''
fastqc_command = ['fastqc', '-o', fastqc_path, '-f',fastq, reads]
subprocess.run(fastqc_command)
    
'''Aligning the reads to the reference genome'''
#cont_outsam = os.path.join(alignment_path, cont_rootname + '.sam')
exp_outsam = os.path.join(alignment_path, exp_rootname + '.sam')

#cont_bowtie_alignments_command = ['bowtie2', '-p', '2', '-x','/Volumes/sesame/joerecovery/genomes/TAIR10_chr_all','-1',cont1,'-2',cont2, '-S',cont_outsam]
exp_bowtie_alignments_command = ['bowtie2', '-p', (int(cores) - 2), '-x',ref,'-1',exp1,'-2',exp2, '-S',exp_outsam]
print('Aligning contuts to reference genome:  ' + str(ref))
# subprocess.run(cont_bowtie_alignments_command)
print('Aligning tests to reference genome:  ' + str(ref))
subprocess.run(exp_bowtie_alignments_command)


'''Continuing the alignment process'''
# cont_outbam = os.path.join(alignment_path, cont_rootname + '.bam')
# cont_bam_command = ['samtools', 'view' ,'-Sb',cont_outsam, '>', cont_outbam]
exp_outbam =  os.path.join(alignment_path, exp_rootname + '.bam')
exp_bam_command = ['samtools', 'view' ,'-Sb',exp_outsam, '>', exp_outbam]
# subprocess.run(cont_bam_command)
subprocess.run(exp_bam_command)

# cont_outsortedbam = os.path.join(alignment_path, cont_rootname + '_sorted.bam')
# cont_sort_command = ['samtools','sort',cont_outbam,cont_outsortedbam]
exp_outsortedbam = os.path.join(alignment_path, exp_rootname + '_sorted.bam')
exp_sort_command = ['samtools','sort',exp_outbam,exp_outsortedbam]
# subprocess.run(cont_sort_command)
subprocess.run(exp_sort_command)

'''Generating the alignment statistics for troubleshooting and bias detection'''
# cont_outflagstat = os.path.join(alignment_path, cont_rootname + '_flagstat.txt')
# cont_flagstat_command = ['samtools', 'flagstat', cont_outsortedbam, cont_outflagstat]
exp_outflagstat = os.path.join(alignment_path, exp_rootname + '_flagstat.txt')
exp_flagstat_command = ['samtools', 'flagstat', exp_outsortedbam, exp_outflagstat]
# subprocess.run(cont_flagstat_command)
subprocess.run(exp_flagstat_command)


'''Reformatting the aligned files to fastq files for transcript quantification'''
# cont_outfastq = os.path.join(alignment_path, cont_rootname + '_bamtofastq.fq')
# cont_bam2fastq_command = ['bedtools', 'bamtofastq', '-i', cont_outsortedbam, '-fq', cont_outfastq]
exp_outfastq = os.path.join(alignment_path, exp_rootname + '_bamtofastq.fq')
exp_bam2fastq_command = ['bedtools', 'bamtofastq', '-i', exp_outsortedbam, '-fq', exp_outfastq]
# subprocess.run(cont_bam2fastq_command)
subprocess.run(exp_bam2fastq_command)




'''Transcript counting using stringtie. This will produce a gtf file of the read
counts for each of the genes involved in the gtf provided.'''
stringtie_out = os.path.join(stringtie_path, exp_rootname + '_read_counts.gtf')
stringtie_command = ["stringtie",  '-eB', '-p', cores, '-G ', stringtie_index, '-o', stringtie_out, exp_outsortedbam]



'''
If there are any major issues with the python scripts, there is a bash script here that can do the same job -
I just find it easier to edit the python scripts because the language is overall much simpler. Modificication of 
the bash script would be needed (altering the paths to the files etc) before this could be run.


set -ueo pipefail

# set the constant variables of your experiment,
# directories for samples, index and the number of threads

IDXGTF=/Users/frankwellmer/refs/Araport11_GFF3_genes_transposons.Mar92021.gff
IDX=~/refs/Hisat2/tair10
SAMPLES=/Volumes/LaCie/READS/LAB_Omics_Data/ALPS_Omics/RNASeq/2021_ChristosVelanis_ALP2-REF6/
THREADS=$1

mkdir -p bam

iter=$(ls ${SAMPLES} | grep fq | sed 's/_[1-2]\.fq.*//g')

for i in ${iter};
do
  hisat2 -p ${THREADS} --max-intronlen 20000 -x ${IDX} -1 ${SAMPLES}${i}_1.fq.gz -2 ${SAMPLES}${i}_2.fq.gz | \
  samtools sort -@ ${THREADS} -o bam/${i}.bam
done

cd bam
ls | parallel samtools index
cd ..


# Estimate transcript abundances and create table counts Ballgown style
for i in ${iter};
do
  stringtie  -eB -p ${THREADS} -G ${IDXGTF} -o ballgown/${i}/${i}.gtf bam/${i}.bam
done

(END)
'''