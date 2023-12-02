import sys
import argparse
import random



'''Small function to randomise a fasta input, often useful to check read quality
after receiving RNA-seq or ChIP-seq sameples'''

'''Setting the arguments - 
Usage = fasta_randomiser.py --fasta FASTAFILE --num Number of sequences you want to take out
'''

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--fasta', required = True, type = str)
parser.add_argument('--num', required =False, default = 50, type = int)
args =parser.parse_args()

fasta = args.fasta
amount = args.num
fasta_outname = fasta + '_random_select'

list = []
r = open(fasta,'r')
newfa = open(fasta_outname,'w+')
fa=r.readlines()
num  = int(len(fa))
rangenum = int(num/amount)

'''Because it is all in a list and we're selecting random numbers there is the
chance that it could be an out of bounds index - just run the script again. Will
try to fix soon'''
for i in range(rangenum-1):
    randnum = random.randint(1, num)
    if randnum % 2 != 0:
        newfa.write(fa[0] + '\n' + fa[randnum] + '\n')
    if randnum % 2 == 0:
        randnum = randnum -1
        newfa.write(fa[0] + '\n' + fa[randnum] + '\n')
