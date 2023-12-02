import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import argparse
import subprocess
seq = 'tcAATCAtac'
stops = ['TAG','TGA','TAA']
promoters = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/PlantPan3.0_promoter_analysis_SPAdes/g29179.t1','r')
seqfile = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/augustus_SPAdes_genomic_assembly/augustu_contigs_Ns.mrna','r')
for seqs in SeqIO.parse(seqfile,'fasta'):
    if seqs.id == 'g18330.t1':
        print(len(seqs.seq))
        gene_with_pro = seqs.seq[491:].upper()
        for i in range(0, len(gene_with_pro)-2, 3):
            codon = gene_with_pro[i:i+3]
            if codon == "ATG":
                start = i
            if codon in stops:
                stop = i
print(gene_with_pro[start:stop])
print(len(gene_with_pro[start:stop]))


print(promoters.head())

tair10 = open('/Users/josephbeegan/Desktop/Work/TAIR10_Chr_all.fa','r')

chrom_loops = pd.read_table('/Users/josephbeegan/Desktop/Work/microarray_SUP/chromatin_loops_weigel.txt')
print(chrom_loops.head())
'''Want to be able to only look at particular chromosomal loops, and then order them based on position'''
chrom_loops_3 = chrom_loops[chrom_loops['Chromosome'] == 3]
chrom_loops_3 = chrom_loops_3.sort_values('Region_a_from')
SUP_postion = 8242256
SUP_end_position = 8243372
chrom_loops_3['Start_distance_from_SUP'] = chrom_loops_3['Region_a_from'] - SUP_postion
chrom_loops_3_5000_from_SUP = chrom_loops_3[chrom_loops_3['Start_distance_from_SUP'].between(-10000,10000)]
print(chrom_loops_3.head(100))
print(chrom_loops_3_5000_from_SUP.head())
print(chrom_loops_3_5000_from_SUP[['Region_a_from','Region_a_to','Region_b_from','Region_b_to']])
print(len(chrom_loops_3_5000_from_SUP))

for i in range(len(chrom_loops_3_5000_from_SUP)):
    print(chrom_loops_3_5000_from_SUP['Region_a_from'].values[i])

for seq in SeqIO.parse(tair10,'fasta'):
    for i in range(len(chrom_loops_3_5000_from_SUP)):
        if str(chrom_loops_3_5000_from_SUP['Chromosome'].values[i]) in seq.id:
            print(seq.id)
            print(seq.seq[chrom_loops_3_5000_from_SUP['Region_a_from'].values[i]:chrom_loops_3_5000_from_SUP['Region_a_to'].values[i]])
            print(seq.seq[chrom_loops_3_5000_from_SUP['Region_b_from'].values[i]:chrom_loops_3_5000_from_SUP['Region_b_to'].values[i]])




outfile = open('/Volumes/LaCie/READS/sinapis_alba/SPAdes_assemblies/full_run_memory_cap_60Gb_merged_not_interlaced/BLAST/global_alignment_spades_assembly_to_tair_mrna_blast.txt','w+')
for seq in SeqIO.parse(contigs,'fasta'):
    for i in range(len(datafram)):
        aligned = pairwise2.align.localxx(seq.seq,datafram['queries'][i])
        for j in aligned:
            outfile.write(seq.id + '\n')
            outfile.write(j)
            outfile.write(datafram['titles'][i])
