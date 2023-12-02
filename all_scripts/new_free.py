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
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import subprocess
import pysam
snapfasta = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_snap_prediction.txt','r')
trinity = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/Trinity.fasta.masked','r')
spades = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdesRNA_assembly.fasta','r')
contigs = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/contigs.fasta','r')
blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/sinapis_augustus_mrna_tairxml','r')
trichome_names_tair = ['AT5G53200', 'AT3G27920','AT1G79840','AT2G46410']
outfile = open('/Users/josephbeegan/Desktop/sinapis_contigs_no_newline.fasta','w+')
soap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_augustus_mrna_blast_TAIR10cDNA.xml','r')
soap_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_augustus.mrna')
spades_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/augustus_SPAdes_genomic_assembly/augustu_contigs_Ns.mrna','r')
spades_snap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_snap_blast_xml','r')
spades_snap_genes = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdes_snap.fasta','r')
soap_snap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_snap_blast_xml','r')

line = ''
bamming_command = ['echo', 'fart', line]
file = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/queries_file.txt','r')
for line2 in file:
    line = line2
    subprocess.Popen(bamming_command, stdout=subprocess.PIPE, bufsize=-1)
