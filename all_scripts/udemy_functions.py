def hello(name):
    print('Howdy ' + name + '!')
    print('Howdy ' + name + '!')


hello('Good')
hello('Gleeb')
hello('I hate you')


def plustOne(number):
    return number + 1


print(plustOne(5))


def div42by(divideBy):
    try:
        return 42 / divideBy
    except ZeroDivisionError:
        print('Error: you tried to divide by zero')
        # This is used to avoid the error popup when dividing
        # by 0
        # Can be useful for input validation


print(div42by(2))
print(div42by(12))
print(div42by(0))
print(div42by(1))

print('How many cats do you have?')
numCats = input()
try:
    if int(numCats) >= 4:
        print('That is a lot of cats')
    elif int(numCats) < 0:
        print('That makes no sense')
    else:
        print('That is not that many cats')
# This is an example of input validation - making sure that the input fits
except ValueError:
    print('You did not enter a number')
