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



file  = '/Volumes/sesame/joerecovery/Project_folder/Nanopore_blast_xmls/TAIR_10_alignment_hits_dedup'
outfile = open('/Volumes/sesame/joerecovery/Project_folder/Nanopore_blast_xmls/TAIR_10_alignment_hits','w+')


count = 0
linecount = 0
with open(file,'r') as read:
    for line in read:
        linecount += 1
        if 'ATC' in line.split(' ')[1]:
            count += 1
        
print(str(count))
print(str(linecount))
print(str(count/linecount*100))
