from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse
import sys
import linecache


my_seq = Seq("TGATCAGT")

def run(args):
    l = Seq(args.input)
    print("Your sequence = " + l)
    print("Reverse = "+l[::-1])
    print("Complement = " + l.complement())

def main():
	parser=argparse.ArgumentParser(description="Trim a reads to a specific length ")
	#Adding arguments is much simpler than expected, the 'dest' component of it
	#is what is called in the functions above, allows for as many arguments input
	#as you  want.
	parser.add_argument("-in",help="Sequence" ,dest="input", type=str, required=True)
	parser.add_argument("-len",help="length of the reads" ,dest="length", type=int, required=False, default=131)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
