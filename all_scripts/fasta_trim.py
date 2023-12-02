import argparse
import linecache
#This is taking an input
#VERY important to have 'run' defined before 'main'
def run(args):
	f = open(args.input)
	l = args.length
	#making the output file (would like to learn how to make this have
	#the same name as the input)
	fout1 = open("trimmed.txt", "w")
	with f as f1:
		with fout1 as f2:
			lines =f1.readlines()
			for line in lines:
				if len(line) <= l:
					f2.write(line + "\n")
				else:
					f2.write(line[0:l] + '\n')
	#Important to close the output file within the defined function

	fout1.close()


def main():
	parser=argparse.ArgumentParser(description="Trim a reads to a specific length ")
	#Adding arguments is much simpler than expected, the 'dest' component of it
	#is what is called in the functions above, allows for as many arguments input
	#as you  want.
	parser.add_argument("-in",help="fasta input file" ,dest="input", type=str, required=True)
	parser.add_argument("-len",help="length of the reads" ,dest="length", type=int, required=False, default=131)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
