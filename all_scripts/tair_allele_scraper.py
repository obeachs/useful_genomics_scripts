from bs4 import BeautifulSoup
import requests
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq

# tair_IDS = []
# tair_CDNA = '/Users/josephbeegan/Desktop/Work/TAIR10_cdna_20101214_updated'
# with open(tair_CDNA,'r') as gene_list:
#     for seq in SeqIO.parse(gene_list,'fasta'):
#         tair_IDS.append(seq.id[0:9])
#
# print(tair_IDS)
#


source = requests.get('https://www.arabidopsis.org/servlet/Search?type=polyallele&action=search&name_type1=locus_name&method1=4&term1=AT1G01010&name_type2=locus_name&method2=2&term2=&allele=any&poly_type=any&insertion_type=any&poly_site=any&inheritance=any&ecoLow=any&ecoHi=any&transgene=any&mutagen=any&size=25&order=name').text
soup = BeautifulSoup(source,'lxml')
results = open('/Users/josephbeegan/Desktop/tair_scraper_results.txt','w+')


table = soup.find_all("table")
table = table[0]
# print(table.prettify)
list = ''
for row in table.find_all('td'):
    list = list + str(row)

#print(list)

if 'deletion' in list:
    print('yay')
else:
    print('farthist')




    headers = []
    for header in table.find_all('th'):
        headers.append(header.text)


    datas = []
    for data in table.find_all('tr'):
        datas.append(data.text)



    new_datas = []
    results.write(final_title + '\n')
    for i in datas:
        new_datas.append(i.split()[1])
    print(new_datas)
    if 'deletion' in new_datas:
        print('Deletion allele present')
    else:
        print('No deletion allele')


    for i in new_datas:
        results.write(i + '\n')

    results.write('\n')

for i in range(10002,1000000):

    source = requests.get('https://www.arabidopsis.org.elib.tcd.ie/servlets/TairObject?id=' + str(i) + '&type=locus').text

    #results = open('/Users/user/Desktop/tair_scraper_results.txt','w+')

    soup = BeautifulSoup(source,'lxml')


    title = soup.find_all('table')
    title = title[0]
    title = title.find('td')
    title = str(title)
    if 'error' in title:
        continue
    else:

        location = title.find('Locus')
        final_title = title[location + 7:location + 16]
        print('Checking TAIR database for : ' + final_title)

        table = soup.find_all("table")
        table = table[13]
        #print(table.prettify)




        headers = []
        for header in table.find_all('th'):
            headers.append(header.text)


        datas = []
        for data in table.find_all('tr'):
            datas.append(data.text)



        new_datas = []
        results.write(final_title + '\n')
        for i in datas:
            new_datas.append(i.split()[1])

        if 'deletion' in new_datas:
            print('Deletion allele present')
        else:
            print('No deletion allele')


        for i in new_datas:
            results.write(i + '\n')

        results.write('\n')
