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


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--input', required = True, default = None, type=str)
args =parser.parse_args()




hic_table = pd.read_table(args.input)
#hic_table = pd.read_table('/Users/josephbeegan/Desktop/Work/microarray_SUP/final_mapped_unique_SRR8089804_1_formatted.fastq_mapped_snippet_1000.txt')
new_columns = ['a','b','c','d','e','f','g','h']
hic_table.columns = new_columns
print(hic_table.a)
print(hic_table.e)
chrom1_1 = hic_table[hic_table.a == 'Chr1']
chrom1_real =chrom1_1[(chrom1_1.a.isin(chrom1_1.e))]
chrom2_1 = hic_table[hic_table.a == 'Chr2']
chrom2_real =chrom2_1[(chrom2_1.a.isin(chrom2_1.e))]
chrom3_1 = hic_table[hic_table.a == 'Chr3']
chrom3_real =chrom3_1[(chrom3_1.a.isin(chrom3_1.e))]
chrom4_1 = hic_table[hic_table.a == 'Chr4']
chrom4_real =chrom4_1[(chrom4_1.a.isin(chrom4_1.e))]
chrom5_1 = hic_table[hic_table.a == 'Chr5']
chrom5_real =chrom5_1[(chrom5_1.a.isin(chrom5_1.e))]



print('Chr1: ' +str(len(chrom1_real.a)))
print('Chr2: ' +str(len(chrom2_real.a)))
print('Chr3: ' +str(len(chrom3_real.a)))
print('Chr4: ' +str(len(chrom4_real.a)))
print('Chr5: ' +str(len(chrom5_real.a)))
