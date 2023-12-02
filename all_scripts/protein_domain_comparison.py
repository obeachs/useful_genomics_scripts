import pandas as pd
import numpy as np
import glob
import boto3


def compare_domain(hitdata_list):
    """Compare the shared and not shared domains in the tables"""
    
    
    
    
protein_domain_data = glob.glob('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/*hitdata*')
for hitdata in protein_domain_data:
    file = pd.read_csv(hitdata, sep='\t',skiprows=7)
    file = file[['Query','Short name','From']]
    grouped_df = file.groupby('Query').apply(lambda x: x.sort_values('From'))

    # reset the index of the resulting dataframe
    file = grouped_df.reset_index(drop=True)
    # create dictionary with gene names as keys and protein domains as values
    gene_domains = {}
    for row in file.itertuples():
        gene_name = row[1][7:]
        domain = row[2]
        gene_domains.setdefault(gene_name, []).append(domain)
    
    # create dictionary with gene names as keys and protein domains as values
    gene_domains_2 = {}
    for hitdata in hitdata_list:
        file = pd.read_csv(hitdata, sep='\t',skiprows=7)
        file = file[['Query','Short name','From']]
        grouped_df = file.groupby('Query').apply(lambda x: x.sort_values('From'))

        # reset the index of the resulting dataframe
        file = grouped_df.reset_index(drop=True)
        # create dictionary with gene names as keys and protein domains as values
        for row in file.itertuples():
            gene_name = row[1][7:]
            domain = row[2]
            gene_domains_2.setdefault(gene_name, []).append(domain)
            
    # create a list of shared domains
    shared_domains = []
    for key in gene_domains:
        if key in gene_domains_2:
            shared_domains.append(gene_domains_2[key])
            for domain in gene_domains[key]:

            
