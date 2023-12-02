import numpy as np
import pandas as pd
import os
import sys
import argparse
import matplotlib.pyplot as plt
import seaborn as sns


#Take input from user to see what samples correspond to each cell
#Samples are in the form of a list of strings
#Config file for this script should be a table delimited file with one column for cell/well 
#coordinates and the other column is the sample names in the format of PRIMER_SAMPLE for experimental 
#samples and REF_SAMPLE for reference samples.


'''Setting up the different variables'''
parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--qpcr', required = True, type = str)
parser.add_argument('--config', required =False, default= 'no', type=str)
args =parser.parse_args()



a = np.array([['a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12'],
['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11','b12'],
['c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12'],
['d1','d2','d3','d4','d5','d6','d7','d8','d9','d10','d11','d12'],
['e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12'],
['f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12'],
['g1','g2','g3','g4','g5','g6','g7','g8','g9','g10','g11','g12'],
['h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12']])

column = []
for i in range(96):
    column.append(i+1)
d = {'col1': ['a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12',
'b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11','b12',
'c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12',
'd1','d2','d3','d4','d5','d6','d7','d8','d9','d10','d11','d12',
'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12',
'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12',
'g1','g2','g3','g4','g5','g6','g7','g8','g9','g10','g11','g12',
'h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12'], 'col2':column}

df = pd.DataFrame(data=d)


qpcr = pd.read_csv(args.qpcr, sep='\t', skiprows=2, header=None)

if '.' in args.qpcr:
    outqpcr = args.qpcr.split('.')[0] + '_qPCR_analysis.tsv'
else:
    outqpcr = args.qpcr + '_qPCR_analysis.tsv'
if '.' in args.qpcr:
    outplt = args.qpcr.split('.')[0] + '_qPCR_analysis.png'
else:
    outplt = args.qpcr + '_qPCR_analysis.png'
qpcr.columns =['TF','colour', 'cell', 'Sample','cp', 'NAME', 'zero','NAME2']
qpcr['cell'] = qpcr['cell'].str.lower()
#Left join df and qpcr dataframes
df = df.merge(qpcr, how='left', left_on='col1', right_on='cell')
#rename 'col1' to 'cell' in df
df = df.drop(columns=['cell'])
#Change NaN to 0 in cp column
df['cp'] = df['cp'].fillna(0)

x = np.array(df['cp'])
shape = (8,12)
print('Preview of plate - cancel script and rearrange coordinates if it looks wrong')
print(x.reshape(shape))



if args.config != 'no':
    file_path = args.config
    if os.path.exists(file_path):
        print('The file exists')
        out_df = pd.read_csv(file_path,sep='\t', header=None)
        print(out_df)
        out_df.columns = ['cell','sam']
        new_df = df.merge(out_df, how='left', left_on='col1', right_on='cell')
        
    else:
        print('The specified file does NOT exist, run the script again with a real config file')
        sys.exit()
    
else:
    number = int(input("Type the number of the sample you want to analyze, including technical replicates and reference samples: "))
    samples = []
    cell_list = []
    for i in range(number):
        sample = input("Enter the name of sample number " + str(i+1) + " to be analyzed. Write in the format of PRIMER_SAMPLE for experimental samples and REF_SAMPLE for reference samples: " )
        samples.append(sample)
        cells = input("Enter cells to be analyzed (separated by commas). For example, if GL1_LEAF was in cells a1 and b1, you type a1,b1 and hit enter: ")
        cell_list.append(cells.split(','))

    out = {'sam':samples, 'cell':cell_list}
    out_df = pd.DataFrame(data=out)

    #Explode the 'cell' column into individual rows
    out_df = out_df.explode('cell')
    new_df = df.merge(out_df, how='left', left_on='col1', right_on='cell')


'''/Users/josephbeegan/qpcr_script_config.txt'''
#Split the 'sam' column into individual columns by separating by '_'
new_df['sam_name'] = new_df['sam'].str.split('_').str[0]
new_df['sam_type'] = new_df['sam'].str.split('_').str[1]

type = new_df['sam_type'].unique().tolist()
type = [x for x in type if str(x) != 'nan']


out_df_list = []
for i in type:
    temp_df = new_df[new_df['sam_type'] == i]
    print('first_temp')
    print(temp_df)
    #temp_df = temp_df[['sam_name','cp']]
    #Group values of temp_df by 'sam_name' and get mean of 'cp' column
    #temp_df = temp_df.groupby('sam_name')['cp'].mean()
    #temp_df = temp_df.groupby(['sam_name']).agg({'cp':['mean','std']})
    if temp_df.empty:
        continue
    refs = temp_df['cp'][temp_df['sam'].str.contains('REF')].mean()
    temp_df = temp_df[~temp_df['sam'].str.contains('REF')]
    temp_df['delta2'] = 2**(refs - temp_df['cp'])
    temp_df = temp_df.groupby(['sam']).agg({'delta2':['mean','std']})
    temp_df = temp_df.reset_index()
    print(temp_df)
    # temp_df = temp_df.reset_index()
    samp = temp_df['sam'].tolist()
    delta2_mean = temp_df['delta2']['mean'].tolist()
    delta2_std = temp_df['delta2']['std'].tolist()
    temp_df_dict = {'sam':samp, 'cp':delta2_mean, 'sd':delta2_std}
    temp_df= pd.DataFrame(data=temp_df_dict)
    temp_df['sam_name'] = temp_df['sam'].str.split('_').str[0]
    temp_df['sam_type'] = temp_df['sam'].str.split('_').str[1]
    out_df_list.append(temp_df)
    # temp_df_dict = {'sam_name':samp, 'cp':cp_mean, 'sd':cp_std}
    # temp_df= pd.DataFrame(data=temp_df_dict)
    # print(temp_df)
    # ref_value = temp_df['cp'][temp_df['sam_name'] == 'REF'].values
    # #Calculate the result column, 2^(ref_value - 'cp')
    # temp_df['result'] = 2**(ref_value - temp_df['cp'])
    # temp_df['sam_type'] = i
    # temp_df = temp_df[temp_df['sam_name'] != 'REF']
    # #Merge 'sam_name' and 'sam_type' columns with '_' as separator
    # temp_df['sample_name'] = temp_df['sam_name'].str.cat(temp_df['sam_type'], sep = '_')
    # #Drop 'sam_type' column
    # temp_df = temp_df.drop(columns=['sam_type','sam_name'])
    # temp_df = temp_df[['sample_name','cp','result','sd']]
    # print(temp_df)
    # out_df_list.append(temp_df)
    # if temp_df.empty:
    #     continue
    # #Print a plot of temp_df with 'cp' on the y-axs using 'std' as the error bars matplotlib
    
    
    #Use seaborn to make a barplot using the temp_df data with error bars for 'result'
    #sns.barplot(x='sample_name', y='result', data=temp_df, color='blue', errcolor='black', yerr=temp_df['sd'])
    #plt.savefig('/Users/josephbeegan/test_out' + str(i) + '.png')   # save the figure to file
    

final_df =  pd.concat(out_df_list)
#Print a plot of final_df with error bars on 'cp' using matplotlib
#yerr=final_df['sd']
# #Use seaborn to make a barplot using the final_df data with error bars for 'result'
sns.set(rc={'figure.figsize':(11.1,11.1)})
#bplot = sns.barplot(x='sam_type', y='cp', data=final_df, color='blue',errcolor='black',errorbar='sd', hue = 'sam_name', palette='magma',capsize=.4,linewidth=3)
#bplot = sns.barplot(data=final_df, x="sam_name", y="cp", hue="sam_type",yerr=final_df['sd'])
#bplot.set_xticklabels(bplot.get_xticklabels(), rotation=90)
#fig = bplot.get_figure()
#fig.savefig(outplt)

# plt.savefig('/Users/josephbeegan/test_out.png')   # save the figure to file


print('Writing file to ' + outqpcr)
final_df.to_csv(outqpcr, sep='\t', index=False, quoting=False)
