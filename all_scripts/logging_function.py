import logging
import pathlib
from pathlib import PurePath
import sys


def sum():
    logging.debug("This is the debug message")
    num1 = 5
    num2 = 10
    sum = num1 + num2
    print(f'Sum:{sum}')

def extractpath(file):
    name = str(pathlib.Path(file).parent.absolute())  + '/'
    return name
def extractitle(file):
    obj = PurePath(file)
    name = obj.name
    return name



def main():


       input1 = sys.argv[1]
       path = extractpath(input1)
       name = extractitle(input1)
       scriptname = extractitle(str(sys.argv[0]))

       fmtstr = "  %(asctime)s: (%(filename)s): %(levelname)s: %(funcName)s  Line: %(lineno)d - %(message)s"
       datestr = "%m/%d/%Y %I:%M:%S %p "

       logging.basicConfig(
       filename=path + name + '_' + scriptname + '_log.txt',
       level=logging.DEBUG,
       filemode="w",
       format=fmtstr,
       datefmt=datestr,
       )

       logging.info("Running " + str(scriptname) + ' on ' + name)

       count = 0
       with open(input1,'r') as fileopen:
           for line in fileopen:
               count += 1

       print(count)



if __name__ == '__main__':
    main()
