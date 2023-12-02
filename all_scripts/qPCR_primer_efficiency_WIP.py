from cmath import log
from tokenize import group
import pandas as pd
import numpy
from scipy.stats import linregress
import argparse

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--table', required = True, type = str)
parser.add_argument('--number', required = False, default = None, type= int)
parser.add_argument('--outfile', required = False, default = None, type = str)
args =parser.parse_args()

tab = args.table

num = 4
if args.number is not None:
    num = args.number

out = tab.split('.')[0] + '_efficiencies_calculated.txt'
if args.outfile is not None:
    out = args.outfile

#Max number of columns for qPCR plate = 12
#Max number of rows for primer efficiencies = 5
#SUM((B5:B9-AVERAGE(B5:B9))*(C5:C9-AVERAGE(C5:C9)))/SUM((B5:B9-AVERAGE(B5:B9))^2)
#Skipping the first row is necessary because of the extra line at the start of 
#Roche files.
tab = pd.read_table(tab, skiprows= [0])

def subset_table(full_table,pairnum):
    potential_pairs = [['1','2','3'],['4','5','6'],['7','8','9'],['10','11','12']]
    subsetted_table = full_table[full_table["Pos"].str[1].isin(potential_pairs[pairnum])]
    return(subsetted_table)

def get_efficiency(subsetted_table):
    grouped = subsetted_table['Cp'].groupby(subsetted_table['Pos'].str[0])
    grouped = list(grouped.mean())
    dilutions =[log(1),log(.2),log(.04),log(.008),log(.0016)]
    slope, intercept, r_value, p_value, std_err = linregress(grouped, dilutions)
    efficiency = ((10**(-1/slope)-1)*100)
    efficiency = str(efficiency).split('+')[0][1:]
    return efficiency

with open(out,'w+') as outfile:
    outfile.write('Primer pair' + '\t' + 'Efficiency' + '\n')
    for i in range(num):
        sub = subset_table(tab,i)
        eff = get_efficiency(sub)
        outfile.write('Primer_' + str(i+1) + '\t' + str(eff) + '%' + '\n')



