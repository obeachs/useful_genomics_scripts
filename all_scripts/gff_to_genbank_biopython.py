import Bio
import sys
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
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
from Bio.Seq import Seq
from BCBio import GFF
from BCBio.GFF import GFFExaminer

def sinapis_id_fix(id):
    start = id.find('Sialb')
    return id[start:start+15]
def tair_id_fix(id):
    start = id.find('AT')
    return id[start:start+9]

gff_file = '/Volumes/sesame/joerecovery/Project_folder//sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Salba_584_v3.1.gene.gff3'
fasta_input = '/Volumes/sesame/joerecovery/Project_folder//sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Phytozome/PhytozomeV13/Salba/v3.1/assembly/Salba_584_v3.0.fa'
out_file = '/Volumes/sesame/joerecovery/Project_folder//sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/Phytozome/PhytozomeV13/Salba/v3.1/assembly/Salba_584_v3.1.gene.gbk'
with open(out_file, "wt") as f_out:
    for record in GFF.parse(gff_file, fasta_input):
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write(record, f_out, "genbank")