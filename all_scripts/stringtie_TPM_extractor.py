import numpy as np
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--gtf', required = False, default = '/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_cdna_20101214_updated', type = str)
args = parser.parse_args()
gtf = args.gtf
outname = gtf.split('.')[0] + '_statlines_only.gtf'


gtf = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/SRR2961888_Salba_new.gtf'
with open(gtf,'r') as file, open(outname,'w+') as out:
    for line in file:
        if 'transcript' in line:
            if 'exon' not in line:
                line = line.replace(";", "\t" )
                out.write(line)