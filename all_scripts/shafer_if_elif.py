language = 'Python'

if language == 'Python':
    print('Condition was True')

user = 'Admin'
logged_in = False

if user == 'Admin' and logged_in:
    print('Admin Page')
else:
    print('Bad')

if not logged_in:
    print('please log')
else:
    print('Welcome')

# Objects, if they have different labels in memory
# a = 1, b = 1 wont be the same object, unless we make
# b=a or vice versa
# To check whether or not they have the same label, use
# print(a is b)
# False Values: False, None, Zero of any type, any empty sequence -(),[],''
# Any empy mapping {}
