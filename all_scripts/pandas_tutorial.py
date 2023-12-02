import numpy as np
from pandas import Series, DataFrame
import pandas as pd

#A normal series in pandas would look like this, can add ,index=['n','n','n']
#to the parentheses to change the index identifying each point.
obj = Series([4,7,-5.3])
print(obj)
obj2 = Series([4, 7, -5, 3], index=['d', 'b', 'a', 'c'])
print(obj2.index)
#Pandas makes it very easy to search for things within the array-like object
#by adding the name of the  index in square parentheses
print(obj2['a'])
#Can apply numpy operations on the array-like objects # TODO:
obj2[obj2>0]
obj2*2
print('b' in obj2)
#Pandas is very good at incorporating data from dictionaries in python
species_data = {"Frogs":4, "Humans":2,"Spiders":6, "Octopus":8}
obj3 = Series(species_data)
#It looks like the name of the object in the dictionary becomes the index of the
#series/object
#Using the Series function also makes it so that the the dictionary is ordered
print(obj3)
#Can be then compared to another list/Series using the indexes
species = ['Humans','Spiders','Slugs','Frogs','Octopus']
#This would be very useful to find out if there is a certain set of genes in
#a large dataset etc.
obj4 = Series(species_data, index=species)
print(obj4)
print(pd.isnull(obj4))
#Dataframes similar to those found in R can be made very easily in pandas
dataset = {"species": species, "legs":[2,6,0,4,8], 'ID':np.arange(5)}
frame = DataFrame(dataset)
print(frame)
print(frame.legs)
#New information can be added to the dataframes (even columns that did not exist beforehand)
frame['eyes'] = [2,8,2,2,2]
print(frame)
#Adding a list to a dataframe will require the list to be the same length
#as the dataframe,  but adding a Series to a dataframe will cause it to conform
#to the dataframe's index
add_to_frame = Series(['Predator','Predator','Predator','Predator'],
index=[0,1,3,4])
frame['PP'] = add_to_frame
print(frame)
print(frame.columns)
