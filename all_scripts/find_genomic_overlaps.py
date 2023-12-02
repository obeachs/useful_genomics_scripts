from email import header
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


sinapis_id = []
tairid = []
tair_start = []
tair_end = []
leng = []
main = []
xml = '/Volumes/sesame/joerecovery/genomes/brapa_R_o_18/Brassica_rapa_ro18.SCU_BraROA_2.3.dna_sm.toplevel.fa_all_ids_BLAST_blast_to_TAIR10_chr_all.fa_xml'
with open(xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        tair_id =alignment.title.split(' ')[0]
                        barley_id = record.query
                        start = hsp.sbjct_start
                        end = hsp.sbjct_end
                        if end > start:
                            length = end-start
                        else:
                            length = start-end
                        sinapis_id.append(barley_id)
                        tairid.append(tair_id)
                        tair_start.append(start)
                        tair_end.append(end)
                        leng.append(length)
                        main.append({'TAIR_ID': tair_id,'Start': start,'End': end,'length':length})

gff = pd.read_table('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff', header=None)
gff = gff.drop(gff[gff[2] =='chromosome'].index)
gff = gff[gff[2] =='gene']
for i in range(len(tairid)):
    if i < 2:
        brapa_start = tair_start[i]
        brapa_end = tair_end[i]
        chr = tairid[i]
        new_gff = gff[gff[0] == chr]
        new_gff['Diff_start_abs'] = abs(new_gff[3] - brapa_start)
        new_gff['Diff_end_abs'] = abs(new_gff[4] - brapa_end)
        new_gff['Diff_start'] = (new_gff[3] - brapa_start)
        new_gff['Diff_end'] = (new_gff[4] - brapa_end)
        new_gff['Length'] = leng[i]
        sorted_start = new_gff.sort_values(by='Diff_start_abs')
        sorted_end = new_gff.sort_values(by='Diff_end_abs')
        final_gff = new_gff
        print(sorted_start.head())
