
from email.quoprimime import quote
from Bio import Restriction as rst
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict
from Bio.Seq import Seq
import re
import pandas as pd
import functools 
import requests
from bs4 import BeautifulSoup




source = requests.get('https://international.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities').text
soup = BeautifulSoup(source,'html.parser')
table = soup.find_all("table")
table = table[0]
    
df = pd.DataFrame(columns=['Enzyme','Restriction_site'])
for row in table.find_all('tr'):    
    # Find all data for each column
    columns = row.find_all('td')
    if(columns != []):
        site = columns[0].text.strip()
        enzyme = columns[1].text.strip()
        df = df.append({'Enzyme':enzyme,'Restriction_site':site}, ignore_index=True)
#df['real_namme'] = df['Enzyme'].str.split(' ')[0]
df['Enzyme'] = df['Enzyme'].apply(lambda x: x.split(' ')[0])
df['Enzyme'] = df['Enzyme'].apply(lambda x: x.split('-')[0])
df['Plain_site'] = df['Restriction_site'].str.replace('/|1|2|3|4|5|6|7|8|9|-','')
df['Plain_site'] = df['Plain_site'].str.replace(r"\(.*\)","")
print(df)

enzymes_list = []
with open('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/KS1-2_genotyping/res_enzymes.txt','r') as restriction:
    for line in restriction:
        name = line.split(' ')[0]
        enzymes_list.append(name)
for i in range(len(enzymes_list)):
    enzymes_list[i] = re.sub('\n','',enzymes_list[i])
print(enzymes_list)
print(df['Enzyme'].to_list())

newdf = df[df['Enzyme'].isin(enzymes_list)]
print(newdf)
df.to_csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/KS1-2_genotyping/res_enzymes_with_sites.csv', index=False)
newdf.to_csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/KS1-2_genotyping/lab_enzymes_with_sites.csv', index=False)



#Function taken from StackOverflow https://stackoverflow.com/questions/20381912/type-object-restrictiontype-has-no-attribute-size/23907326#23907326

def create_enzyme(name):
    e_types = [x for t, (x, y) in typedict.items() if name in y][0]
    enzyme_types = tuple(getattr(rst, x) for x in e_types)
    return rst.RestrictionType(name, enzyme_types, rest_dict[name])

seq = Seq('CCTATTTCCAACCCTTTCTCCTCCATCCTCACCAAG')
a = rst.Analysis(rst.AllEnzymes, seq)
#rb = functools.reduce(lambda x, y: x + y, map(create_enzyme, enzymes_list))

alist = a.full() 
