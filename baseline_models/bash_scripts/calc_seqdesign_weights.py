#!/usr/bin/env python
import os
import glob
import numpy as np
from tqdm.contrib import tmap
from collections import defaultdict

alphabet = 'ACDEFGHIKLMNPQRSTVWY'
theta = 0.2

os.makedirs("../../alignments/wavenet/DMS/", exist_ok=True)

for a2m_filename in glob.glob("../../alignments/focus_column_only/DMS/*.a2m"):
    a3m_filename = a2m_filename.replace(".a2m", ".a3m").replace("focus_column_only", "full")
    output_filename = a2m_filename.replace(".a2m", ".fa").replace("focus_column_only", "wavenet")
    print(a2m_filename)

    aa_dict = {}
    for i,aa in enumerate(alphabet):
        aa_dict[aa] = i

    seq_name_to_sequence = defaultdict(str)
    seq_names = []

    name = ''
    INPUT = open(a2m_filename, 'r')
    counter = 0
    for i, line in enumerate(INPUT):
        line = line.rstrip()
        if line.startswith('>'):
            name = line.split(' ')[0]
            seq_names.append(name)
            #print name
            counter += 1
        else:
            seq_name_to_sequence[name] += line
    INPUT.close()

    focus_seq_name = seq_names[0]

    # Select focus columns
    focus_seq = seq_name_to_sequence[focus_seq_name]
    focus_cols = [ix for ix, s in enumerate(focus_seq) if s == s.upper()]
    focus_seq_trimmed = [focus_seq[ix] for ix in focus_cols]
    seq_len = len(focus_cols)
    alphabet_size = len(alphabet)


    seq_name_to_focus_col_sequence = {}
    for seq_name,sequence in seq_name_to_sequence.items():
        # Replace periods with dashes (the uppercase equivalent)
        sequence = sequence.replace('.','-')

        #then get only the focus columns
        seq_name_to_focus_col_sequence[seq_name] = [sequence[ix].upper() for ix in focus_cols]

    # Remove sequences that have bad characters
    alphabet_set = set(list(alphabet))
    seq_names_to_remove = []
    for seq_name,sequence in seq_name_to_focus_col_sequence.items():
        for letter in sequence:
            if letter not in alphabet_set and letter != '-':
                seq_names_to_remove.append(seq_name)

    seq_names_to_remove = list(set(seq_names_to_remove))

    for seq_name in seq_names_to_remove:
        del seq_name_to_focus_col_sequence[seq_name]

    # Encode the sequences
    print ('Encoding sequences', flush=True)
    x_train = np.zeros((len(seq_name_to_focus_col_sequence.keys()),len(focus_cols),len(alphabet)))
    x_train_name_list = []
    for i,seq_name in enumerate(seq_name_to_focus_col_sequence.keys()):
        sequence = seq_name_to_focus_col_sequence[seq_name]
        x_train_name_list.append(seq_name)
        for j,letter in enumerate(sequence):
            if letter in aa_dict:
                k = aa_dict[letter]
                x_train[i,j,k] = 1.0
    print("sequences encoded", flush=True)

    list_seq = x_train
    list_seq = list_seq.reshape((list_seq.shape[0], list_seq.shape[1] * list_seq.shape[2]))
    def compute_weight(seq):
        number_non_empty_positions = np.dot(seq,seq)
        if number_non_empty_positions>0:
            denom = np.dot(list_seq,seq) / np.dot(seq,seq)
            denom = np.sum(denom > 1 - theta)
            return 1/denom
        else:
            return 0.0 #return 0 weight if sequence is fully empty
    weights = np.array(list(tmap(compute_weight, list_seq, total=list_seq.shape[0])))

    print("Neff:",np.sum(weights))

    name_to_weight = {x_train_name_list[i]:weights[i] for i in range(len(x_train_name_list))}


    def seq_newline_formatter(seq,num_char=60):
        seq_format = ''
        for i in range(0,len(seq),num_char):
            seq_format += seq[i:i+num_char]+'\n'
        return seq_format

    lc_sequence_name_to_sequence = {}
    INPUT = open(a3m_filename, 'r')
    seq = ''
    name = ''
    first_time = True
    for line in INPUT:
        line = line.rstrip()
        if line[0] == '>' and first_time:
            first_time = False
            name = line.split(' ')[0]
        elif line[0] == '>' and first_time == False:
            lc_sequence_name_to_sequence[name] = seq
            name = line.split(' ')[0]
            seq = ''
        else:
            seq_chunk = ''
            for aa in line:
                if aa != '.' and aa != '-':
                    seq_chunk += aa.upper()
            seq += seq_chunk
    lc_sequence_name_to_sequence[name] = seq
    INPUT.close()

    counter = 0
    with open(output_filename, 'w') as OUTPUT:
        for name,weight in name_to_weight.items():
            if name in lc_sequence_name_to_sequence:
                OUTPUT.write(name+":"+str(weight)+'\n'+seq_newline_formatter(lc_sequence_name_to_sequence[name]))
                counter += 1
            else:
                print("skipped", name)

    print("Num sequences written out:",counter, flush=True)
