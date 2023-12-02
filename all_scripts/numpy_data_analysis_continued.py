import numpy as np
#Creating a new array by taking the value from X if
#its corresponding value in cond is TRUE, otherwise
#taking it from Y - the following is the slower pure
#python method. DOES NOT WORK IN MULTIDIMENSIONAL
#ARRAYS.
xarray = np.array([1.1,1.2,1.3,1.4,2.5])
yarray = np.array([2.1,2.2,2.3,2.4,2.5])
cond = np.array([True, False, True, True, False])
r = [(x if c else y) for x,y,c in zip(xarray, yarray, cond)]
print(r)
#Much easier and more efficient to do this with np.where
r2 = np.where(cond, xarray, yarray) #the second and third
#arguments here dont need to be arrays.
print(r2)
#np.where is generally used to create a new array based on
#another array.
wherearray = np.random.randn(4,4)
print(wherearray)
#Here is a method to make all positive values with 2 and all
#negative values with -2
print(np.where(wherearray > 0,2, -2))
#this will change only the positive values to 2. We can think
#of np.where as an ifelse kind of thing: np.where(if,do x, else do y)
print(np.where(wherearray > 0, 2, wherearray))
#Having .mean at the end of the arrayname, will calculate the mean
#of all of the elements in the arrays.
wherearray2 = np.where(wherearray > 0,2,3 )
print(wherearray2)
print(wherearray2.mean())
#Can specify what row or axis to do the statistical function on
print(wherearray2.mean(axis=1))
print(wherearray2.sum(0))
#Sort function could come in handy for data analysis - can sort
#the whole array by value or sort each invdividual dimension of values
#by value by passing an axis number with sort
sortarray = np.random.randn(5,3)
print(sortarray)
print(np.unique(sortarray))
print(sortarray.sort(1))
