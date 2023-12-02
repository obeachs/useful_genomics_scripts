import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from Bio import Align
parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--queryfasta', required = True, default = None, type = str, help = 'Homer-annotated ChIP-seq peak file')
args =parser.parse_args()
aligner = Align.PairwiseAligner()



query = args.queryfasta
out = query.split('.')[0] + '_fixed.fa'
with open(query,'r') as que, open(out, 'w+') as outfile:
    seen = []
    count = 1
    for seq in SeqIO.parse(que,'fasta'):
        count += 1
        if seq.id not in seen:
            seen.append(seq.id)
            outfile.write('>' + str(seq.id) + '\n' + str(seq.seq)+ '\n')
        else:
            new_id = seq.id + '_' + str(count)
            outfile.write('>' + str(new_id) + '\n' + str(seq.seq)+ '\n')

#         print('Aligning  ' + str(seq.seq) + '\n' + 'to  ' + str(main_seq))
#         alignments = aligner.align(main_seq,str(seq.seq))
#         alignments = alignments[0]
#         score = alignments.score
#         scores.append(score)

# df_list = [ids,scores]
# df = pd.DataFrame(list(zip(ids, scores)), columns=['ID', 'Score'])
# df = df.sort_values('Score')
# print(df)


