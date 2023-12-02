import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2
import random


parser.add_argument('--sequence', required = False, default = 'All', type = str)
args =parser.parse_args()


if args.sequence == 'All'
    x= input("Enter variables: ").split(",")
