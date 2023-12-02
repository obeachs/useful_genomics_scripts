import numpy as np
from pandas import Series, DataFrame
import pandas as pd
#May also be useful to create or encounter data that is in the form of a nested
#dictionary
pop = {'Nevada': {2001: 2.4, 2002: 2.9},'Ohio': {2000: 1.5, 2001: 1.7, 2002: 3.6}}
popframe = DataFrame(pop)
print(popframe)
#If we want to transpose this dataframe - make the index the columns and vice
#versa, we simply do: popframe.T
print(DataFrame(pop, index = [2000, 2001, 2004]))
#Possibly important to understand the nature of arithmetic methods in
#pandas
df1 = DataFrame(np.arange(15).reshape((3,5)), columns = list('abcde'))
df2 = DataFrame(np.arange(12).reshape((3,4)), columns = list('abcd'))
#If we simply do df2 +df1 then we will get NaNs because they are different
#sizes - which could be what we want since it is only adding the cells
#that contain values
#We need to do df2.add(df1, fill_value=0)
print(df2.add(df1, fill_value=0))
SUP = pd.read_csv("~/Desktop/Work/microarray_SUP/temp_SUP/Final_data/merge3.5D_sorted_fix_29-01-20_0.05_test.txt")
print(SUP)
SUP_wanted = SUP[["locus","logFC","adj.P.Val"]]
print(SUP_wanted)
del SUP['t']
del SUP['Annotation']
print(SUP.columns)
print(SUP[SUP_wanted['logFC'] < -1])
