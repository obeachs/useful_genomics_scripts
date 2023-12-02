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
import requests
import itertools
import argparse
import subprocess
from os import listdir
from os.path import isfile, join
import glob

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--query', required = True, default = None, type = str)
parser.add_argument('--blastdatabase', required = True, default = None, type= str)
parser.add_argument('--outfolder', required = True, default = None, type = str)
args =parser.parse_args()
query = args.query
blastdatabase= args.blastdatabase

'''Just a function to check the first line of a file'''
def read_first_lines(filename, limit):
  result = []
  with open(filename, 'r') as input_file:
    # files are iterable, you can have a for-loop over a file.
    for line_number, line in enumerate(input_file):
      if line_number > limit:  # line_number starts at 0.
        break
      result.append(line)
  return result

print(query)



if str(args.outfolder)[-1] != '/':
    args.outfolder = str(args.outfolder + '/')
queryname = query.split('/')[-1]
queryname = queryname.split('.')[0]
outname = args.outfolder + queryname + '_blast_to_' + blastdatabase.split('/')[-1] +  '_xml'
print(outname)
blast_command = ['blastn', '-query', query, '-db',  blastdatabase, '-out', outname, '-max_target_seqs','6', '-num_threads', '8', '-evalue', '1e-6', '-outfmt', '5']
subprocess.run(blast_command)






# '''Using the read_first_lines function to make sure that the database files were are using are fasta files'''
# full_path_file_list = []
# file_list = []
# test_count = 3
# count = 0
# for f in listdir(blastdb_folder):
#     print(f)
#     print(blastdb_folder+'/'+f)
#     filename = blastdb_folder+'/'+f
#     if 'n' not in filename[-3:]:
#         firstline = read_first_lines(filename,1)
#         if ">" in firstline[0]:
#             file_list.append(f)
#             full_path_file_list.append(blastdb_folder+'/'+f)
# print(file_list)
    
