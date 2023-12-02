import sys
import Bio
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import UnknownSeq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten
from tensorflow.keras.models import Sequential
from sklearn.metrics import confusion_matrix
import itertools
import tensorflow.keras.backend as K
import argparse
import scipy


def fasta_sequences_recover(fasta):
    sequence_list = []
    for seq in SeqIO.parse(fasta,'fasta'):
        seqs = str(seq.seq)
        sequence_list.append(seqs)
    return(sequence_list)

def sequences_to_one_hot(sequence_array):
  input_sequences = []
  integer_encoder = LabelEncoder()
  one_hot_encoder = OneHotEncoder(categories = 'auto')
  for sequences in sequence_array:
    integer_encoded = integer_encoder.fit_transform(list(sequences))
    integer_encoded = np.array(integer_encoded.reshape(-1,1))
    one_hot_encoded = one_hot_encoder.fit_transform(integer_encoded)
    input_sequences.append(one_hot_encoded.toarray())
  return(input_sequences)


def labels_to_one_hot(inputs):
  one_hot_encoder = OneHotEncoder(categories='auto')
  labels = np.array(inputs).reshape(-1,1)
  input_labels = one_hot_encoder.fit_transform(labels).toarray()
  return(input_labels)
