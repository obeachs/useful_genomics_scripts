import argparse
import linecache
def run(arg):
	qual = arg.quality_score
	file = open(arg.input) #The names of the in and out files will be given when the script is called
	output = open(read_1.fastq)
	read_1.fastq = []
	with open(file) as f:
		for i, line in enumerate(f,1):
			if "1:N:0" in linecache.getLine(file, i-1):
				linestart = i
				lineend = i + 3
				for n in range(linestart, lineend):
					read_1.append(linecache.getline(file, j))
	read_1.close()	
def(main):	
	parser=argparse.ArgumentParser(description="Split a paired-end fastq file into 2 seperate read files")
	parser.add_argument("-in", help="Paired-end fastq input file", dest="input")
	parser.set_defaults(func=run)
	arg=parser.parse_args()
	arg.func(args)

if __name__=="__main__":
	main()