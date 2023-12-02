import numpy as np
names = np.array(['Bob','Joe','Will','Bob','Will','Joe','Joe'])
data = np.random.randn(7,4)
print(data)
#We can compare the array to a string to get a
#boolean array
print(names == 'Bob')
#Using the resulting boolean array, we can
#index another array's values. Important to
#note that the indexed array must be the same
#length as the boolean array - creates a copy
print(data[names != "Bob"])
#Can also manipulate arrays by comparing to
#boolean array.
data[names != "Will"] = 7
print(data)
#'Fancy indexing'
for i in range(7):
    data[i] = i
#We can reshape the arange array into whatever
#shape we want, from 1D to multidimensional
array = np.arange(32).reshape((8,4))
#This will now be indexed to select very specific
#element from the array with one line of code.
#Taking (1,0),(5,3),(7,1) and (2,2) all in one go
#Each list of numbers is a different axis
print(array[[1,5,7,2],[0,3,1,2]])
#Keep in mind that fancy indexing, unlike slicing,
#always copies the data into a new array.
#If using a list to select rows from an array, negative
#indices will select from the end
print("TRANSPOSING AND SWAPPING \n")
#Adding .T at the end of an array will transpose it.
#Transposing is making columns rows and rows columns
newarray = np.arange(15).reshape((3,5))
print(newarray)
print(newarray.T)
#Transposing with multidimensional arrays is more complicated
#All of the below are views of the data made, not overwriting or
#copying it.
multiarray = np.arange(16).reshape((2,2,4))
print(multiarray.transpose((1,0,2)))
print(multiarray.swapaxes(1,2))
