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


nanopore_reads = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Nanopore_first_run_6-9-21_sinapis_alba.fasta'

tair_xml = '/Volumes/sesame/joerecovery/Project_folder/Nanopore_blast_xmls/TAIR10_cdna_20101214_updated_Nanopore_first_run_6-9-21_sinapis_alba_xml'

chloroplast_xml = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/chloroplast_assembly_analysis/Nanopore_first_run_to_chloroplast_xml'

table = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/chloroplast_assembly_analysis/reads_and_starting_points_dedup')
table.sort_values(by=['24875'])
print(table['24875'])
tair_chloroplast_hits = []
with open(tair_xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if 'C' in alignment.title.split(' ')[1]:
                            tair_chloroplast_hits.append(record.query)

start_list = []
with open(chloroplast_xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        start_list.append(hsp.sbjct_start)
                    
                        
start_list.sort()
print(start_list)