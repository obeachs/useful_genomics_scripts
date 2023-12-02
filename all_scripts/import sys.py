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

blastx = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/swissprot_blastx/all_rna_seq_merged_maingenome.fa_all_ids_BLAST_blast_to_uniprot_swissport_full.fasta_xml'

def blast_xml_analysis(xml):
    # out.write('sinapis_id' + '\t' + 'tair_id' + '\t' + 'sinapis_start' + '\t' + 'sinapis_end' + '\t' + 'hit_length' + '\t' + 'tair_start' + '\t'+'tair_end' + '\n')
    with open(xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        print("DESCRIPTION")
                        print(description.title)
                        print('ALIGNMENT')
                        print(alignment.title)
blast_xml_analysis(blastx)