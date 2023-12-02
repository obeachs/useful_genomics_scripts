import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2
import glob




parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--IDS', required = False, default = 'All', type = str)
parser.add_argument('--cdna',required = True, default = None, type = str)
parser.add_argument('--db', required = True, default = None, type= str)
parser.add_argument('--out', required = True, default = None, type = str)
parser.add_argument('--blast', required = False, default = 'n', type = str)
args =parser.parse_args()
cdna = args.cdna
out = args.out
blastdb = args.db
blast = args.blast
if out[-1] != '/':
    out = out + '/'
print(out)
subprocess.run(['mkdir','-p', out])

blast_command = ['python3.10', '/Volumes/sesame/joerecovery/scripts/full_blast_analysis_script_full.py','--cdnafile',cdna, '--blastdatabase', blastdb, '--outfolder', out, '--blast', blast]
subprocess.run(blast_command)

os.chdir(out)
blast_output = ''
for file in glob.glob('*dedup.txt'):
    blast_output = file

parsing_command = ['python3.10', '/Volumes/sesame/joerecovery/scripts/parsing_blast_analysis.py', '--blastoutput',blast_output,'--queryfasta', cdna]
subprocess.run(parsing_command)
subprocess.run(['mkdir', '-p', 'alignments'])

parsed_output = ''
for file in glob.glob('*parsed.fa'):
    parsed_output = out + file

print(parsed_output)
os.chdir('alignments')

align_command = ['python3.10','/Volumes/sesame/joerecovery/scripts/generate_fasta_for_alignments.py','--fasta', parsed_output]
subprocess.run(align_command)
