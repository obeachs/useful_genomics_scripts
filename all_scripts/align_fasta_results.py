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
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import Align
from Bio import pairwise2



fasta_1 = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/blast_script_runs/SOAP_denovo_assembly/contigs_with_several_hits_dedup.fasta'
fasta_2 = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/blast_script_runs/SPAdes_ragtag_with_read_correction_BO/contigs_with_several_hits_dedup.fasta'




for seq in SeqIO.parse(fasta_1,'fasta'):
    for seq2 in SeqIO.parse(fasta_2,'fasta'):
        aligner = Align.PairwiseAligner()
        outfile = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/blast_script_runs/seq.id + '_' + seq2.id + "_alignment.txt"
        alignment = aligner.align(seq.seq,seq2.seq)
