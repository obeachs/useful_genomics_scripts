import numpy as np
from pandas import Series, DataFrame
import pandas as pd
KNU5 = pd.read_csv("~/Desktop/Work/microarray_KNU_timecourse/results/KNU_timecourse_5D_all_sorted.txt")
KNU5.name = "KNU5"
KNU6 = pd.read_csv("~/Desktop/Work/microarray_KNU_timecourse/results/KNU_timecourse_6D_all_sorted.txt")
KNU6.name = "KNU6"
KNU7 = pd.read_csv("~/Desktop/Work/microarray_KNU_timecourse/results/KNU_timecourse_7D_all_sorted.txt")
KNU7.name = "KNU7"
print(KNU5.columns)
for frame in (KNU5,KNU6,KNU7):
    frame_new_cols = frame
    new_col = ["GeneName",frame.name+"logFc",frame.name+"AveExpr",frame.name+"adj.P.Val",frame.name+"ext_gene",frame.name+"Annotation" ]
    frame_new_cols.columns = new_col
    print(frame.columns)
#It seems that you can only merge one thing at a time
#It also doesn't like when you merge things with the same
#column names - causes some issues.
#This is a method which combines the two datasets by changing
#the column names bar one, which could be useful in some ways
#but probably not very useful otherwise
print(KNU5.columns)
shared_5D_6D = KNU5.merge(KNU6)
print(shared_5D_6D)
print(shared_5D_6D.drop_duplicates(subset="GeneName"))
shared_5D_6D.to_csv(r'~/Desktop/PANDAS_5D6D_different_columns.csv')

KNU5d = pd.read_csv("~/Desktop/Work/microarray_KNU_timecourse/results/KNU_timecourse_5D_all_sorted.txt")
KNU6d = pd.read_csv("~/Desktop/Work/microarray_KNU_timecourse/results/KNU_timecourse_6D_all_sorted.txt")
KNU7d = pd.read_csv("~/Desktop/Work/microarray_KNU_timecourse/results/KNU_timecourse_7D_all_sorted.txt")
shared_5D_6D_same_columns = KNU5d.append(KNU6d, ignore_index= True)
shared_5D_6D_same_columns_noNA = shared_5D_6D_same_columns.dropna()
print(shared_5D_6D_same_columns.count())
print(shared_5D_6D_same_columns_noNA.count())
#Reason for this being so long is that the append function doesn't remove duplicates
#so we have to remove them manually
##Lines below are just ways to look at the correct length of the data array
print(KNU5d.shape[0])
print(shared_5D_6D_same_columns_noNA.shape[0])
shared_5D_6D_same_columns_noNA.to_csv(r'~/Desktop/PANDAS_5D6D')
shared_5D_6D_same_columns_noNA = shared_5D_6D_same_columns_noNA.drop_duplicates()
print(shared_5D_6D_same_columns_noNA.shape[0])
#Can print more than one column at a time
print(KNU5[["KNU5logFc","KNU5ext_gene"]])
print(KNU5.KNU5logFc)
#This will print the first values for all columns in the dataframe
print(KNU5.loc[0])
#Easy way to filter out certain data
top_5D = (KNU5['KNU5logFc'] > 2)
print(KNU5[top_5D])
#We cannot use AND and OR in PANDAS, we use & and | respectively
top_5D_modal = (KNU5['KNU5logFc'] > 2) | (KNU5['KNU5logFc'] < -2)
print(KNU5.loc[top_5D_modal, ["GeneName",'KNU5ext_gene','KNU5logFc']])
#This might be a useful way to find overlaps between datasets without having to
#merge them and then trim later
sevendaygenes = KNU7['GeneName']
shared_5_7 = KNU5["GeneName"].isin(sevendaygenes)
print(KNU5.loc[shared_5_7, "GeneName"])
