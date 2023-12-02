import sys
import argparse


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--inputfile', required = True, type=str)
parser.add_argument('--outputfile', required = True, type=str)
args = parser.parse_args()
input_file = open(args.inputfile,'r')
output_file = args.outputfile
lines_seen = set() # holds lines already seen
with open(output_file, "w+") as output_file:
	for each_line in input_file:
	    if each_line not in lines_seen: # check if line is not duplicate
	        output_file.write(each_line)
	        lines_seen.add(each_line)
