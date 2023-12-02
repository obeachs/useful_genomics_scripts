import sys
import re
import itertools
import argparse
import subprocess
import os


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--infile', required = True, default = None, type=str)
args =parser.parse_args()

file = args.infile
print(str(file))
new_file_name = str(file) + '_clean_text'
print(new_file_name)


def isLineEmpty(line):
    return len(line.strip()) == 0


with open(file) as fileoopen:
    lines = fileoopen.readlines()
    for i in range(len(lines)):
        if lines[i].startswith('Query'):
            start_table = i
            break
    with open(new_file_name,'w+') as newfilewrite:
        for line in lines[start_table:]:
            if ' - >' in line:
                line = line.replace(' - >','\t')
                newfilewrite.write(line)
            else:
                newfilewrite.write(line)

with open(new_file_name,'r') as filehandle:
        lines = filehandle.readlines()

with open(new_file_name, 'r') as filehandle:
    newer_filename = new_file_name + '_finished'
    with open(newer_filename,'w+') as newer_file2:
        lines = filter(lambda x: x.strip(), lines)
        for i, line in enumerate(lines):
            newer_file2.write('%d%s'%(i, line))
