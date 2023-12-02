import random
from bs4 import BeautifulSoup
import requests
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq


tair_IDS = []
tair_CDNA = '/Users/frankwellmer/joe/conorSS/sequences/TAIR10_seq_20101214_updated'
with open(tair_CDNA,'r') as gene_list:
    for seq in SeqIO.parse(gene_list,'fasta'):
        tair_IDS.append(seq.id[0:9])
print(tair_IDS[0:14])


results_list_current = []
results = '/Users/frankwellmer/Desktop/tair_scraper_results.txt'
with open(results,'r') as res:
	for line in res:
		results_list_current.append(line.split(':')[0])

print(results_list_current[0:14])

def non_match_elements(list_a, list_b):
    non_match = []
    for i in list_a:
        if i not in list_b:
            non_match.append(i)
    return non_match

new_tair_IDS = non_match_elements(tair_IDS, results_list_current)
print(new_tair_IDS[0:14])


number_list = []
number = len(new_tair_IDS)
for x in range(number):
	w = random.randint(1,number)
	if w not in number_list:
		number_list.append(w)
	else:
		continue

with open(results,'a') as res:
	for x in number_list:
		name = new_tair_IDS[x]
		source = requests.get('https://www.arabidopsis.org/servlet/Search?type=polyallele&action=search&name_type1=locus_name&method1=4&term1=' + name + '&name_type2=locus_name&method2=2&term2=&allele=any&poly_type=any&insertion_type=any&poly_site=any&inheritance=any&ecoLow=any&ecoHi=any&transgene=any&mutagen=any&size=25&order=name').text
		soup = BeautifulSoup(source,'lxml')
		table = soup.find_all("table")
		table = table[0]
		# print(table.prettify)
		list = ''
		for row in table.find_all('td'):
			list = list + str(row)

		if 'deletion' in list:
			print('yay')
			print(name)
			res.write(str(name) + ': Deletion allele found' + '\n')
		else:
			continue
