import numpy as np
import pandas as pd

def rearrange_hitdata(table):

    outputname = table.split('/')[-1][:-4]
    outputname = outputname + '_rearranged.txt'
    outputpath = table2.split('/')[:-1]
    outputpath = '/'.join(outputpath)
    outputpath = outputpath + '/' + outputname
    table = pd.read_table(table)
    table = table[['Query','Short name']]
    print(outputpath)

    list_for_empty_df = ['Query','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25']
    empty = pd.DataFrame(columns=['Query', '', 'Action'])

    setname = ''
    list_of_ID_names = []
    for index, row in table.iterrows():
        if row['Query'] != setname:
            setname = str(row['Query'])
            list_of_ID_names.append(setname)

    ID_list = table['Query'].tolist()

    max_number = 0
    for i in list_of_ID_names:
        if ID_list.count(i) > max_number:
            max_number= ID_list.count(i)
    list_for_empty_df_new = list_for_empty_df[:max_number +1]
    print(list_for_empty_df_new)
    emptydf = pd.DataFrame(columns=list_for_empty_df_new)
    print(emptydf)


    biglist = []
    for title in list_of_ID_names:
        list = []
        for index, row in table.iterrows():
            if row['Query'] == title:
                list.append(row['Short name'])
        list.insert(0,title)
        print(str(len(list)))
        biglist.append(list)


    df = pd.DataFrame(biglist,columns = list_for_empty_df_new)
    df.to_csv(outputpath, sep = '\t',index = False)


'''Very easy script to use, just type the full file path in the function below,
it will create a file called FILE_rearranged.txt'''
rearrange_hitdata('/Volumes/seagatedrive/downloads/TRY_hitdata(1).txt')
