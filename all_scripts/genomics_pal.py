from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd





parser = argparse.ArgumentParser(description = 'ChipSeq prep for inputs and for test files')
parser.add_argument('--file1', required = True, default = None, type = str)
parser.add_argument('--file2', required = True, default = None, type= str)
parser.add_argument('--input1', required = True, default = None, type= str)
parser.add_argument('--input2', required = True, default = None, type= str)
parser.add_argument('--fastqc', required = False, default = 1, type= str)
parser.add_argument('--hisat2', required = True, default = 0, type= str)
parser.add_argument('--bowtie', required = True, default = 1, type= str)
parser.add_argument('--filter', required = True, default = 1, type= str)
parser.add_argument('--outdir', required = True, default = 1, type= str)

args =parser.parse_args()


def title_maker(string1, string2):
    if '/' in string1:
        string1 = string1.split('/')[-1]
        

list_of_files


def fastqc_run(reads):  
    fatqc_command = ['fastqc' reads '-o' args.outdir]


def align_reads_bowtie(read1,read2):
    align_command = ['bowtie2' '-p' 6 '-q' '--local' '-x' $genome '-1' read1 '-2' read2 '-S'   ]




blast_command = ['blastn', '-query', query, '-db',  db, '-out', outname, '-max_target_seqs','6', '-num_threads', '8', '-evalue', '1e-6', '-outfmt', '5']
    subprocess.run(blast_command)