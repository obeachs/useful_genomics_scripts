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

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--queryIDS', required = False, default = 'All', type = str)
args =parser.parse_args()

def blast_xml_analysis(xml):
    with open(xml,'r') as xmlblast:
        soaptair_record = NCBIXML.parse(xmlblast)
        descriptions = []
        for record in soaptair_record:
            for description in record.descriptions:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        tair_id =alignment.title.split(' ')[0]
                        barley_id = record.query
                        print(barley_id)
                        if record.query == barley_id:
                            print(tair_id)
                        else:
                            barley_id = record.query
                            print(barley_id)
                            print(tair_id)
                        continue
blast_xml_analysis(xml)