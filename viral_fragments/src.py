import sys
import os
import numpy as np 
import pandas as pd

import plotly.express as px
import plotly.graph_objects as go

UKN = -1000

from sklearn.metrics import roc_auc_score

import torch

'''
    Functions for fragmenting RNA sequences based on DMS data

    >>> fragment_RNA('AAAAAAAUUUUUUU', [[0,5],[1,4], [8,10]], [0,0,1,1,0,0,1,1,0,1,0,1,1,1], 'test', min_length=3, min_unpaired_length=1, min_auroc=0.8)

'''

def fragment_RNA(sequence, paired_bases, signal, name, data_type, min_length=100, min_unpaired_length=10, min_auroc=0.8):

    paired_bases = np.array(paired_bases)
    signal = np.array(signal)

    assert len(signal) == len(sequence), 'Sequence and signal must be the same length'

    # Find all non paired bases to get candidate cutting points
    all_paired_bases = np.sort(paired_bases.flatten())

    len_nonPaired_regions = np.diff(all_paired_bases) - 1
    cut_points = np.round( (all_paired_bases[:-1] + all_paired_bases[1:]) / 2).astype(int)

    # Find unstructured regions
    dot = np.array(['.']*len(signal))
    dot[paired_bases[:,0]] = '('
    dot[paired_bases[:,1]] = ')'
    dot = ''.join(dot)

    unstructured_regions = np.zeros(len(signal), dtype=bool)
    bracket_counter = 0
    for i in range(len(dot)):

        if dot[i] == '(':
            bracket_counter += 1
        elif dot[i] == ')':
            bracket_counter -= 1
        else:
            if bracket_counter == 0:
                unstructured_regions[i] = True

        assert bracket_counter >= 0

    # Convert to binary for AUROC
    isUnpaired = np.ones(len(signal))
    isUnpaired[paired_bases.flatten()] = 0


    # Find cutting points
    
    cut_idxs = []

    idx_start = -1
    for i in range(len(cut_points)):

        if (unstructured_regions[cut_points[i]] 
            and (cut_points[i] - idx_start > min_length) 
            and len_nonPaired_regions[i] > min_unpaired_length):

            dms_window = signal[idx_start+1:cut_points[i]]
            pair_window = isUnpaired[idx_start+1:cut_points[i]][dms_window>=0]
            if len(np.unique(pair_window))==2:
                auroc = roc_auc_score(pair_window, dms_window[dms_window>=0])
                # print(auroc)

                if auroc>min_auroc:
                    cut_idxs.append(cut_points[i])
                    idx_start = cut_points[i]
    
    # Last segment
    if len(signal) - idx_start > min_length:
        dms_window = signal[cut_points[-1]+1:]
        pair_window = isUnpaired[cut_points[-1]+1:][dms_window>=0]

        if len(np.unique(pair_window))==2:
            auroc = roc_auc_score(pair_window, dms_window[dms_window>=0])

            if auroc>min_auroc:
                cut_idxs.append(len(signal))

    ## Output dataframe
    data_struct = {}
    start_idx = -1
    for i, end_idx in enumerate(cut_idxs):

        sub_seq = sequence[start_idx+1:end_idx+1]
        sub_dms = signal[start_idx+1:end_idx+1]
        sub_struct = dot[start_idx+1:end_idx+1]

        data_struct[name+'_'+str(i)] = {'sequence': sub_seq, data_type: sub_dms.tolist(), 'structure': sub_struct}
        start_idx = end_idx

    return pd.DataFrame(data_struct).T
    # import json 
    # json.dump(data_struct, open('data/sars_dms_fragments.json', 'w'), indent=2)



def compare_frags(sequences, signals, paired_bases, min_f1=0.9, name='', min_length=100, min_unpaired_length=10, min_auroc=0.8):

    max_len = min([len(seq) for seq in sequences])
    
    matrix_1 = pairList2pairMatrix(paired_bases[0], max_len)
    matrix_2 = pairList2pairMatrix(paired_bases[1], max_len)

    idx_to_keep = [[] for _ in range(len(sequences))]
    frag_dfs = []

    for i, (seq, signal, paired_base) in enumerate(zip(sequences, signals, paired_bases)):
        
        frag_df = fragment_RNA(seq, paired_base, signal, 
                               name=name, min_length=min_length, min_unpaired_length=min_unpaired_length, min_auroc=min_auroc)
        frag_dfs.append(frag_df)

        for idx, row in frag_df.iterrows():

            sub_seq = row['sequence']
            i_start, i_end = seq.find(sub_seq), seq.find(sub_seq) + len(sub_seq)

            f1 = compute_f1( matrix_1[i_start:i_end, i_start:i_end], matrix_2[i_start:i_end, i_start:i_end] )
            if f1 > min_f1:
                idx_to_keep[i].append(idx)


    best_i = np.argmax([len(idx_to_keep[i]) for i in range(len(sequences))])
    print('Best fragment: ', best_i)
    return frag_dfs[best_i].loc[idx_to_keep[best_i]]


def compute_f1(pred_matrix, target_matrix, threshold=0.5):
    """
    Compute the F1 score of the predictions.

    :param pred_matrix: Predicted pairing matrix probability  (L,L)
    :param target_matrix: True binary pairing matrix (L,L)
    :return: F1 score for this RNA structure
    """

    pred_matrix = (pred_matrix > threshold).float()

    sum_pair = torch.sum(pred_matrix) + torch.sum(target_matrix)

    if sum_pair == 0:
        return 1.0
    else:
        return (2 * torch.sum(pred_matrix * target_matrix) / sum_pair).item()
    
def pairList2pairMatrix(pair_list, len_seq):
    pair_list = np.array(pair_list).astype(int)
    pairing_matrix = torch.zeros((len_seq, len_seq), dtype=torch.bool)

    if len(pair_list) > 0:
        pairing_matrix[pair_list[:,0], pair_list[:,1]] = 1
        pairing_matrix[pair_list[:,1], pair_list[:,0]] = 1

    return pairing_matrix