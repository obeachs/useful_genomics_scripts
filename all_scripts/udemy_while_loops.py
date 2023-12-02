spam = 0
while spam < 5:
    print('hello world')
    spam = spam + 1


print('My name is')
for i in range(5):
    print('Jimmy FIVE TIMES ' + str(i))
for i in range(5, -1, -1):
    print('JIMMY FIVE TIMES ' + str(i))
for i in range(0, 10, 2):
    print('JIMMY FIVE TIMES ' + str(i))

# Gauss problem, adding all the numbers from 1-100
total = 0
for num in range(101):
    total = total + num
print(total)

import random

print(random.randint(1, 1421049109401))

import pyperclip
