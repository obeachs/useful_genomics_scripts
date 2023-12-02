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
import matplotlib.pyplot as plt
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2




parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--fasta', required = False, default = '/Volumes/sesame/joerecovery/genomes/TAIR10_cdna_20101214_updated', type = str)
args = parser.parse_args()

outfile = args.fasta.split('.')[0] + '_tair_sequences_removed.fa'


with open(args.fasta,'r') as fa, open(outfile,'w+') as out:
    for seq in SeqIO.parse(fa,'fasta'):
        if len(seq.id) != 9 and seq.id[0:1] !='AT':
            out.write('>' + str(seq.id) + '\n' + str(seq.seq) +'\n')