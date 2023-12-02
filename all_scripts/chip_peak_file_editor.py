import sys
import Bio
import re
import os
import numpy as np
import pandas as pd
import requests
import itertools
import argparse
import subprocess
import multiprocessing


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--peakfile', required = True, default = None, type = str)
args =parser.parse_args()


file = args.peakfile

peakfile = pd.read_table(file)

outfile = file.split('.')[0]
outfile = outfile + '_homer_edited.broadPeak'

'''Naming and reordering the columns to match the homer annotation format'''
peakfile.columns = ['Chromosome','Start','End','Peak_ID','Score','Strand','Signal_value','p-value','q-value']

cols = peakfile.columns.tolist()

newlistorder = [3,0,1,2,5,4,6,7,8]

newcols = [cols[i] for i in newlistorder]


peakfile = peakfile[newcols]
peakfile['Chromosome'] = peakfile['Chromosome'].str.lower()

'''Writing the new file out to a new table based on the given filename'''
peakfile.to_csv(outfile,sep = '\t',index = False)
