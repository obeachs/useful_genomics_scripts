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
import itertools
import argparse
import subprocess

table = pd.read_table('/Volumes/seagatedrive/downloads/hitdata_full_bnapus_clean2.txt')


print(table.columns)
print(table['Query'][1])

list_of_genes = ['CDX79408','CDY53885']

#Here is the function to condense the tables into the different genes
#You should pass through the dataframes you made in R as the table variable, this will be with all of the genes from all species
#Then, you need to make lists of genes from each species - i.e AGAMOUS_genes = ['tair_gene_ID','brassica_gene_ID',etc,etc]
#Once you have the lists and the full data from you can pass them through this function
# Example: AGAMOUS_protein_domain_table = create_condensed_data_frames(table, AGAMOUS_genes)
#Then you can output the new table using different methods that are easily found online
def create_condensed_data_frames(dataframe,list):
    new_df = pd.DataFrame(columns=['Query', 'Hit type', 'PSSM-ID', 'From', 'To', 'E-Value', 'Bitscore', 'Accession', 'Short name', 'Incomplete', 'Superfamily'])
    for gene in list:
        title = gene + '_table'
        title = dataframe[dataframe['Query'].str.contains(gene)]
        new_df = new_df.append(title, ignore_index=True)
    return(new_df)

new_table = create_condensed_data_frames(table,list_of_genes)
print(len(table))
print(len(new_table))

new_table.to_csv('/Users/josephbeegan/Desktop/testy_table.txt')
