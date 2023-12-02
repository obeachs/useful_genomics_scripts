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




transcript = '/Volumes/sesame/joerecovery/sinapis_assembly_shenanigans/nanopore_reads_corrected/gene_prediction/RNAseq_to_nanopore_reads_5000_blast_to_TAIR10_cds_20101214_updated_dedup.txt'
genomics = '/Volumes/sesame/joerecovery/sinapis_assembly_shenanigans/nanopore_reads_corrected/all_nanopores_combind5000_k17k25k31.fa'

