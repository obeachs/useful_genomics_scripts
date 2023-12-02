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
from macpath import sep
from pickle import FALSE



ngd = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/slim.ngd')
gff = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/TAIR_gff.gff')





ngd["START"] = ""
ngd['STOP'] = ""
ngd['CHR'] = ''

print(ngd.head())
print(gff.head())

for index, row in ngd.iterrows():
    for fart, poop in gff.iterrows():
        if row['V4'] == poop['V4']:
            row['START'] = poop['start']
            row['STOP'] = poop['end']
            row['CHR'] = poop['chr'].lower()
            print(row['V4'])
      
            
print(ngd.head())
print(gff.head())


ngd.to_csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/fixed_ngd.ngd',sep = '\t',index = False,quoting = None)



