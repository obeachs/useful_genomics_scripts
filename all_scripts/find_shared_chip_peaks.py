import sys
import os
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
parser.add_argument('--peakfileone', required = True, default = None, type=argparse.FileType('r'))
parser.add_argument('--peakfiletwo', required = True, default = None, type=argparse.FileType('r'))
parser.add_argument('--outdir', required = False, default = None, type=str)
args =parser.parse_args()
def listToString(s):

    # initialize an empty string
    str1 = ""

    # traverse in the string
    for ele in s:
        str1 += ele

    # return string
    return str1
def change_gene_name(column_list):
    for n, i in enumerate(column_list):
        if i == 'Gene Name':
            column_list[n]='Gene_Name'
            return(column_list)

one_set = args.peakfileone
two_set = args.peakfiletwo
title_one_set = str(one_set.split('/')[-1])
title_one_set_no_ext = str(title_one_set.split('.')[0])
cwd = os.getcwd()
if args.outdir is not None:
    outdir = args.outdir
else:
    outdir = cwd
title_two_set = str(two_set.split('/')[-1])
title_two_set_no_ext = str(title_two_set.split('.')[0])
UNO = pd.read_table(two_set)
DOS = pd.read_table(one_set)
#for column in UNO.columns:
#    column = list(column)
#    for n, i in enumerate(column):
#        if i == ' ':
#            column[n]='_'
#            column2 = listToString(column)
#            list2.append(column2)
#        else:
#            column2 = listToString(column)
#            list2.append(column2)
#UNO.columns = list2
print(UNO.columns)

column_namesUNO = list(UNO.columns)
column_namesDOS = list(DOS.columns)
#print(column_names)
#for n, i in enumerate(column_names):
#    if i == 'Gene Name':
#        column_names[n]='Gene_Name'
#print(column_names)
column_namesUNO = change_gene_name(column_namesUNO)
column_namesDOS = change_gene_name(column_namesDOS)

UNO.columns = column_namesUNO
DOS.columns = column_namesDOS
print(UNO.columns)
print(DOS.columns)


shared_genes = list(set(UNO['Gene_Name'].unique()) & set(DOS['Gene_Name'].unique()))
print(shared_genes)
shared_genes = [x for x in shared_genes if str(x) != 'nan']
print(shared_genes)
print(list(UNO['Gene_Name']))

print("Found " + str(len(shared_genes)) + " gene peaks between " + title_one_set + ' and ' + title_two_set)
UNOnew = UNO[UNO['Gene_Name'].isin(shared_genes)]
DOSnew = DOS[DOS['Gene_Name'].isin(shared_genes)]
print(len(UNOnew['Gene_Name']))
print(len(DOSnew['Gene_Name']))
print('Found ' + str(len(shared_genes)) + ' shared_genes')
dupesUNO = UNOnew[UNOnew.duplicated(['Gene_Name'])]
print('One-peaks written to ' + outdir+title_two_set_no_ext+ '_shared_with_' + title_one_set_no_ext + '.txt')
dupesUNO.to_csv(outdir+title_two_set_no_ext+'_shared_with_'+title_one_set_no_ext+'.txt',sep='\t')
dupesDOS = DOSnew[DOSnew.duplicated(['Gene_Name'])]
print('Two-peaks written to ' + outdir+title_one_set_no_ext+'_shared_with_'+title_two_set_no_ext+'.txt')
dupesDOS.to_csv(outdir+title_one_set_no_ext+'_shared_with_'+title_two_set_no_ext+'.txt',sep = '\t')
