import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import r2_score

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from seismic2dreem.dump import dump_json
import tqdm

from datetime import datetime

UKN = -1000

dirname = os.path.dirname(__file__)


def import_CT(path_to_ct):

    valid_bases = ['A', 'C', 'G', 'U']

    # Read all CT files to get the list of sequences and pairing matrices
    data = {}
    for file in os.listdir(path_to_ct):
        if file.endswith(".txt") or file.endswith(".ct"):
            with open(os.path.join(path_to_ct, file)) as f:
                data_file = f.readlines()

                sequence = []
                base_pair = []
                for row in data_file[1:]:

                    if row=='':
                            continue
                    
                    row_data = row.split()
                    
                    if row_data[1] == 'T' or row_data[1] == 't':
                        sequence.append('U') # replace T by U
                    else:
                        sequence.append(row_data[1].upper())
                    if int(row_data[4]) != 0:
                        base_pair.append(sorted([int(row_data[0])-1, int(row_data[4])-1])) # Sort and -1 to get 0-based indexing

                # Check if the sequences and structures are valid
                base_pair = np.unique(base_pair, axis=0)
                if len(base_pair) > 0:
                    if (all(b in valid_bases for b in sequence) 
                        and (base_pair<len(sequence)).all() 
                        and len(np.unique(base_pair[:,0]))==len(base_pair)
                        and len(np.unique(base_pair[:,1]))==len(base_pair) ):
                        data[file.split('.')[0]] = {'sequence': ''.join(sequence), 'paired_bases': base_pair.tolist()}
                else:
                    if (all(b in valid_bases for b in sequence)):
                        data[file.split('.')[0]] = {'sequence': ''.join(sequence), 'paired_bases': []}


    # Postprocess data to remove duplicate sequences with different structures
    df = pd.DataFrame.from_dict(data, orient='index')

    unique_seqs = {}
    for i, seq in enumerate(df.sequence):

        struct = np.array(df['paired_bases'][i])

        if seq not in unique_seqs.keys():
            unique_seqs[seq] = [struct]
        else:
            unique_seqs[seq].append(struct)

    # Only keep sequence duplicate that have the same pairing matrix
    for seq in unique_seqs.keys():
        if len(unique_seqs[seq]) > 1:
            matches = []
            for i in range(len(unique_seqs[seq])):
                for j in range(i+1, len(unique_seqs[seq])):
                    matches.append(np.array_equal(unique_seqs[seq][i], unique_seqs[seq][j]))

            # Remove sequences from dataframe if f1 score is not 1.0, keep only one of the structures otherwise
            if not np.array(matches).all():
                df = df[df['sequence'] != seq]
            else:
                row = df[df['sequence'] == seq].head(1)
                df = df[df['sequence'] != seq]
                df = pd.concat([df, row])

    return df


def get_ref2plate(path_excel: str):

    '''
    Get a dictionary that maps reference name to plate number

    Parameters:
        - path_excel (str): Path to the excel file containing the list of references in order

    Returns:
        - ref2plate (dict): Dictionary that maps reference name to plate number

    Example:
    >>> path_excel = 'data/2021-03-03_oligos.xlsx'
    >>> ref2plate = get_ref2plate('test_files/S029472-R1_SynbioTechnologies.xlsx')
    >>> print(ref2plate)
    {'ENSG00000000005.6': 0, 'ENSG00000001036.14': 0, 'ENSG00000001561.7': 0, 'ENSG00000001629.10': 0}
    
    '''
    # Read all sequences
    ref_excel = []
    for ref in pd.read_excel(path_excel)['Oligo Name']:
        if isinstance(ref, str):
            if '_fwd' in ref:
                ref_excel.append(ref.split('_fwd')[0])

    ref2plate = {}
    for i, ref in enumerate(ref_excel):

        plate96 = i//96

        if plate96 in [32, 33, 46, 47]:
            ref2plate[ref] = 3

        elif plate96 in range(12, 16):
            ref2plate[ref] = 11
        
        elif plate96 < 32:
            ref2plate[ref] = plate96//4

        else:
            ref2plate[ref] = (plate96-2)//4

        # print(i, ref, plate96, ref2plate[ref], sep='\t')

    return ref2plate


def read_seismic_output(output_path: str, min_cov=1000, ref2plate: dict=None):

    """
    Read seismic output and return a dataframe with all the information
    Filter bases with low coverage

    Parameters:
        - output_path (str): Path to the output folder of seismic. 
        Should be arranged as output_path/table/sample_name/reference_name/full/mask-per-pos.csv
        with sample_name starting with 'Un' for untreated and '1-5' for treated

        - ref2plate (dict): Dictionary that maps reference name to plate number
        - min_cov (int): Minimum number of bases covered for a base to be kept

    Returns:
        - seismic_df (pd.DataFrame): Dataframe with all the information from seismic output

    Example:
    >>> output_path = 'test_files/test_samples'
    >>> ref2plate = {'ref1': 0, 'ref2': 1}
    >>> seismic_df = read_seismic_output(output_path, min_cov=1000, ref2plate=ref2plate)
    Number of references in output:  3
    >>> true_df = pd.DataFrame({'sample': ['1-5_1_A', '1-5_1_B--Un_2', '1-5_1_B--Un_2'], \
                                'reference': ['ref1', 'ref2', 'ref1'], \
                                'plate': [1, 2, 1], \
                                'replicate': ['A', 'Untreated', 'B'], \
                                'sequence': ['CCTA', 'AUUC', 'CCTA'], \
                                'sub_rate': [[0.0125 , 0.0125, 0.0125, UKN],        [0.0025, 0.0025, 0.0025, 0.025],    [0.0125 , 0.0125, 0.0125, UKN]], \
                                'coverage': [[1000 , 1000, 1000, UKN],              [1000, 1000, 1000, 1000],           [1000 , 1000, 1000, UKN]], \
                                'sub_A': [[0.0025252525252525255, 0, 0.0125, UKN],              [0, 0, 0.0025, 0.012658228],        [0.0025252525252525255, 0, 0.0125, UKN]], \
                                'sub_C': [[0, 0, 0, UKN],                             [0.0025, 0, 0, 0],                  [0, 0, 0, UKN]], \
                                'sub_G': [[0.0025252525252525255, 0.0125, 0, UKN],        [0, 0.0025, 0, 0.012658228],        [0.0025252525252525255, 0.0125, 0, UKN]], \
                                'sub_T': [[0.005037783375314861, 0, 0, UKN],        [0, 0, 0, 0],                       [0.005037783375314861, 0, 0, UKN]], \
                                'n_reads': [14, 14, 14] })
    >>> pd.testing.assert_frame_equal(seismic_df.sort_values(by=['reference', 'replicate']).reset_index(drop=True), \
                                        true_df.sort_values(by=['reference', 'replicate']).reset_index(drop=True)) 

    """


    ## Import from Seismic directly

    # Parse output file to find mask-per-pose.csv files in output_path/table
    #seismic_df = pd.DataFrame(columns=['sample', 'reference', 'plate', 'replicate', 'sequence', 'sub_rate', 'coverage', 'sub_A', 'sub_C', 'sub_G', 'sub_T', 'n_reads'])
    lines = []
    n_references = 0

    # Go over all seismic output
    table_path = os.path.join(output_path, 'table')
    for sample_name in tqdm.tqdm(os.listdir(table_path), desc='Samples', total=len(os.listdir(table_path))):
        sample_path = os.path.join(table_path, sample_name)
        if not os.path.isdir(sample_path):
            continue

        # Extract plate and replicate from sample name
        n_plates = []
        replicates = []
        for s in sample_name.split('--'):
            if s.startswith('Un'):
                n_plates.append(int(s.split('Un_')[1].split('_')[0])-1)
                replicates.append("Untreated")
            elif s.startswith('1-5'):
                n_plates.append(int(s.split('1-5_')[1].split('_')[0])-1)
                replicates.append(s.split('1-5_')[1].split('_')[1])

        # Go through all references       
        for ref_name in tqdm.tqdm(os.listdir(sample_path), desc="Reading through sample {}".format(sample_name), total=len(os.listdir(sample_path))):
            ref_path = os.path.join(sample_path, ref_name, 'full')
            if not os.path.isdir(ref_path):
                continue
            
            n_references += 1
            
            if not ref2plate is None:
                # Reference should be in the right plate
                if not ref2plate[ref_name] in n_plates:
                    continue

            for root, dirs, files in os.walk(ref_path):
                for file in files:
                    if file == 'mask-per-pos.csv':
                        
                        # Extract coverage and sub_rate
                        df = pd.read_csv(os.path.join(root, file), index_col=0)
                        # df = df.loc[:, ['Base', 'Covered', 'Matched', 'Mutated']]

                        coverage = df['Covered'].to_numpy()/4

                        sub_rate = {}
                        for sub_type in ['Mutated', 'Subbed-A', 'Subbed-C', 'Subbed-G', 'Subbed-T']:
                            sub_rate[sub_type] = (df[sub_type]/(df[sub_type]+df['Matched'])).to_numpy()
                            sub_rate[sub_type][np.isnan(sub_rate[sub_type])] = UKN
                            sub_rate[sub_type][coverage<min_cov] = UKN

                        coverage[np.isnan(coverage)] = UKN
                        coverage[coverage<min_cov] = UKN

                        n_reads = len(pd.read_csv(os.path.join(root, 'mask-per-read.csv.gz'), compression='gzip'))/4
                        replicate = replicates[0] if ref2plate is None else replicates[np.where(ref2plate[ref_name]==np.array(n_plates))[0][0]]
                        plate = n_plates[0] if ref2plate is None else ref2plate[ref_name]+1

                        # Save to seismic_df
                        # lines.append([ sample_name, ref_name, 
                        #                                     plate, replicate,
                        #                                     ''.join(df['Base'].tolist()), 
                        #                                     sub_rate['Mutated'].tolist(), 
                        #                                     coverage.tolist(), 
                        #                                     sub_rate['Subbed-A'].tolist(),
                        #                                     sub_rate['Subbed-C'].tolist(),
                        #                                     sub_rate['Subbed-G'].tolist(),
                        #                                     sub_rate['Subbed-T'].tolist(),
                        #                                     n_reads])
                        lines.append({
                            'sample': sample_name,
                            'reference': ref_name,
                            'plate': plate,
                            'replicate': replicate,
                            'sequence': ''.join(df['Base'].tolist()),
                            'sub_rate': sub_rate['Mutated'].tolist(),
                            'coverage': coverage.tolist(),
                            'sub_A': sub_rate['Subbed-A'].tolist(),
                            'sub_C': sub_rate['Subbed-C'].tolist(),
                            'sub_G': sub_rate['Subbed-G'].tolist(),
                            'sub_T': sub_rate['Subbed-T'].tolist(),
                            'n_reads': n_reads
                        })
                        
    seismic_df = pd.DataFrame(lines)

    # Duplicate checks
    # assert seismic_df.duplicated(subset=['reference']).sum() == seismic_df.duplicated(subset=['sequence']).sum(), 'There are duplicated sequences or references'
    # assert seismic_df.duplicated(subset=['reference', 'replicate']).sum() == 0, 'There are duplicated references and replicates'  
    # assert np.all([rep in ['A', 'B', 'Untreated'] for rep in seismic_df['replicate'].unique()]), 'There are replicates other than A, B and Untreated'  
    print('Number of references in output: ', n_references)

    return seismic_df

def filter_bad_references(seismic_df, min_reads:int=1000, signal_tresholds:tuple=(0.3, 0.7), min_covered:float=0.5, primers=None):

    '''

    Filter out references that don't meet certain criteria

    Parameters:
        - seismic_df (pd.DataFrame): Dataframe with all the information from seismic output
        - min_reads (int): Minimum number of reads for a reference to be kept
        - min_cov (int): Minimum number of bases covered for a reference to be counted
        - min_covered (float): Minimum fraction of bases covered for a reference to be kept

    Returns:
        - filtered_df (pd.DataFrame): Filtered dataframe with all the information from seismic output

    Example:
    # ref1 has high signal in Untreated (not consistent mutation) -> all ref1 are removed
    # ref2 has B with too few reads, Untreated with high signal (update sequence) -> A and Untreated are kept, with new sequence
    # ref3 is removed because Untreated has signal in between bounds
    # ref4 is removed because of not well covered
    >>> seismic_df = pd.DataFrame({'sample': ['1_A', '1_B', '1_U', '2_A', '2_B', '2_U', '2_B', '2_U', '3_A' ], \
                                'reference': ['ref1', 'ref1', 'ref1', 'ref2', 'ref2', 'ref2', 'ref3', 'ref3', 'ref4'], \
                                'plate': [1, 1, 1, 2, 2, 2, 2, 2, 3], \
                                'replicate': ['A', 'B', 'Untreated', 'A', 'B', 'Untreated', 'B', 'Untreated', 'Untreated'], \
                                'sequence': ['UAGCU', 'UAGCU', 'UAGCU', 'UCGCU', 'UCGCU', 'UCGCU', 'UAUAU', 'UAUAU', 'UGAAUCU'], \
                                'sub_rate': [[0.0, 0.1 , 0.01, 0.1, 0.0], [0.0, 0.1 , 0.01, 0.1, 0.0], [0.0, 0.8 , 0.01, 0.01, 0.0],        [0.0, 0.1 , 0.01, 0.1, 0.0], [0.0, 0.1 , 0.01, 0.1, 0.0], [0.0, 0.01 , 0.01, 0.75, 0.0],        [0.0, 0.1 , 0.01, 0.1, 0.0], [0.0, 0.01 , 0.01, 0.5, 0.0],          [0.0, 0.01 , UKN, UKN, UKN, 0.1, 0.0]], \
                                'coverage': [[1000, 1000, 1000, 1000], [1000, 1000, 1000, 1000, 1000], [1000, 1000, 1000, 1000, 1000],      [1000, 1000, 100, 1000, 1000], [1000, 1000, 1000, 1000, 1000], [1000, 1000, 1000, 1000, 1000],  [1000, 1000, 1000, 1000, 1000], [1000, 1000, 1000, 1000, 1000],     [1000, 1000 , 1000, 1000, 1000, 1000, 1000]], \
                                'sub_A':    [[None, None, None, None], [None, None , None, None, None], [0.0, 0.0 , UKN, 0.01, 0.0],        [None, None, None, None, None], [None, None, None, None, None], [0, 0, UKN, 0.75, 0],           [None, None, None, None, None], [0, 0, UKN, 0.5, 0],                [UKN, UKN, UKN, UKN, UKN, 0.1, 0]], \
                                'sub_C':    [[None, None, None, None], [None, None , None, None, None], [0.0, 0.3 , UKN, 0, 0.0],           [None, None, None, None, None], [None, None, None, None, None], [0, 0, UKN, 0, 0],              [None, None, None, None, None], [0, 0, UKN, 0, 0],                  [UKN, UKN, UKN, UKN, UKN, 0, 0]], \
                                'sub_G':    [[None, None, None, None], [None, None , None, None, None], [0.0, 0.2 , UKN, 0, 0.0],           [None, None, None, None, None], [None, None, None, None, None], [0, 0, UKN, 0, 0],              [None, None, None, None, None], [0, 0, UKN, 0, 0],                  [UKN, UKN, UKN, UKN, UKN, 0, 0]], \
                                'sub_T':    [[None, None, None, None], [None, None , None, None, None], [0.0, 0.3 , UKN, 0, 0.0],           [None, None, None, None, None], [None, None, None, None, None], [0, 0, UKN, 0, 0],              [None, None, None, None, None], [0, 0, UKN, 0, 0],                  [UKN, UKN, UKN, UKN, UKN, 0, 0]], \
                                'n_reads': [100, 1000, 1000, 1000, 100, 1000, 1000, 1000, 1000], \
                                 })


    >>> filtered_df = filter_bad_references(seismic_df, min_reads=1000, signal_tresholds=(0.3, 0.7), min_covered=0.5, primers=('U','U'))
    -> Before filtering: Number of sequences: 9 and number of unique sequences: 4
    -> After filtering low reads: Number of sequences: 7 and number of unique sequences: 4
    -> After filtering high mutations: Number of sequences: 3 and number of unique sequences: 2
    (1 references were updated and 2 references were removed)
    -> After filtering not well covered: Number of sequences: 2 and number of unique sequences: 1

    >>> true_df = pd.DataFrame({'sample': ['2_A', '2_U'], \
                                'reference': ['ref2', 'ref2'], \
                                'plate': [2, 2], \
                                'replicate': ['A', 'Untreated'], \
                                'sequence': ['UCGAU', 'UCGAU'], \
                                'sub_rate': [[UKN, 0.1 , UKN, 0.1, UKN], [UKN, 0.01 , UKN, 0, UKN]], \
                                'coverage': [[UKN, 1000, UKN, 1000, UKN], [UKN, 1000, UKN, 1000, UKN]], \
                                'n_reads': [1000, 1000] })
    
    >>> pd.testing.assert_frame_equal(filtered_df, true_df) 

    '''
    assert signal_tresholds[0] < signal_tresholds[1], 'signal_thresholds[0] should be smaller than signal_thresholds[1]'
    filtered_df = seismic_df.copy()

    # Filter out low reads references
    print(f'-> Before filtering: Number of sequences: {len(filtered_df)} and number of unique sequences: {len(filtered_df.sequence.unique())}')
    filtered_df = filtered_df[filtered_df['n_reads']>=min_reads].reset_index(drop=True)
    print(f'-> After filtering low reads: Number of sequences: {len(filtered_df)} and number of unique sequences: {len(filtered_df.sequence.unique())}')

    # Filter out bases of the primers
    if primers:
        assert len(primers)==2, 'Primers should be a tuple of 2 strings'
        assert (len(primers[0]) != 0) and (len(primers[1]) != 0), 'Primers should not be empty'

        for i, row in filtered_df.iterrows():
            assert row['sequence'].startswith(primers[0]) and row['sequence'].endswith(primers[1]), 'Sequence does not start or end with the primers'

            columns = ['sub_rate', 'coverage', 'sub_A', 'sub_C', 'sub_G', 'sub_T']
            for col in columns:
                array = np.array(row[col])
                array[:len(primers[0])] = UKN
                array[-len(primers[1]):] = UKN
                filtered_df.at[i, col] = array.tolist()
        

    # Filter out references with high signals in untreated
    ref_to_remove = []
    n_updated = 0
    for ref in filtered_df[filtered_df['replicate']== 'Untreated']['reference'].unique():
        ref_df = filtered_df[(filtered_df['reference']==ref) & (filtered_df['replicate']=='Untreated')]

        sub_rate = np.array(ref_df['sub_rate'].tolist()[0])
        if ( (sub_rate[sub_rate!=UKN] >= signal_tresholds[0]) & (sub_rate[sub_rate!=UKN] <= signal_tresholds[1]) ).any():
            ref_to_remove.append(ref)
            continue

        mask_high_mut = sub_rate > signal_tresholds[1]
        if mask_high_mut.any():

            all_sub_types = np.vstack([ref_df[sub_type].tolist()[0] for sub_type in ['sub_A', 'sub_C', 'sub_G', 'sub_T']]).T
            sorted_all_sub_type = np.sort(all_sub_types, axis=1)

            # If all the mutation are always the same, update the sequence
            if (sorted_all_sub_type[mask_high_mut][:, -1] > 10*sorted_all_sub_type[mask_high_mut][:, -2]).all():
                sequence = np.array([*ref_df['sequence'].tolist()[0]])
                dms = np.array(ref_df['sub_rate'].tolist()[0])
                sequence[mask_high_mut] = np.array(['A', 'C', 'G', 'T'])[np.argmax(all_sub_types[mask_high_mut], axis=1)]
                dms[mask_high_mut] = 0

                assert filtered_df.loc[filtered_df['reference']==ref, 'sequence'].iloc[0] != ''.join(sequence)
                filtered_df.loc[filtered_df['reference']==ref, 'sequence'] = ''.join(sequence)
                filtered_df.at[filtered_df.index[(filtered_df['reference']==ref) & (filtered_df['replicate']=='Untreated')][0], 'sub_rate'] = dms.tolist()
                n_updated += 1

            # Otherwise there might be several versions of the gene, remove the reference
            else:
                ref_to_remove.append(ref)

    filtered_df = filtered_df[~filtered_df['reference'].isin(ref_to_remove)]
    print(f'-> After filtering high mutations: Number of sequences: {len(filtered_df)} and number of unique sequences: {len(filtered_df.sequence.unique())}')
    print(f'({n_updated} references were updated and {len(ref_to_remove)} references were removed)')
    
    # Filter out references that are not well covered
    low_covered = []
    for i, row in filtered_df.iterrows():
        seq = np.array([*row['sequence']])
        mask_AC = (seq=='A') | (seq=='C')
        dms = np.array(row['sub_rate'])
        cov = np.array(row['coverage'])

        # At least a fraction of min_covered of A and C should be covered
        if sum(dms[mask_AC] != UKN)/mask_AC.sum() < min_covered: 
            low_covered.append(i)   
        
        # Filter out G/C bases
        dms[~mask_AC] = UKN
        cov[~mask_AC] = UKN
        filtered_df.at[i, 'sub_rate'] = dms.tolist()
        filtered_df.at[i, 'coverage'] = cov.tolist()

    filtered_df = filtered_df.drop(low_covered)     
    print(f'-> After filtering not well covered: Number of sequences: {len(filtered_df)} and number of unique sequences: {len(filtered_df.sequence.unique())}')

    return filtered_df.drop(columns=['sub_A', 'sub_C', 'sub_G', 'sub_T'], axis=1).reset_index(drop=True)


def get_plate_info(seismic_df):
    
    for plate in seismic_df['plate'].unique():
        plate_df = seismic_df[seismic_df['plate']==plate]

        print("\nPlate", plate)
        print(  "Number of A replicates: ", len(plate_df[plate_df['replicate']=='A']), 
                " | Number of B replicates: ", len(plate_df[plate_df['replicate']=='B']), 
                " | Number of untreated replicates: ", len(plate_df[plate_df['replicate']=='Untreated']) )

        
        # Get all DMS signal of untreated data
        all_dms_untreated = np.concatenate(plate_df[plate_df['replicate']=='Untreated']['relate_sub_rate'].values)
        all_dms_untreated = all_dms_untreated[all_dms_untreated!=UKN]
        all_dms_treated = np.concatenate(plate_df[plate_df['replicate']!='Untreated']['relate_sub_rate'].values)
        all_dms_treated = all_dms_treated[all_dms_treated!=UKN]

        print("Min treated DMS: ", np.min(all_dms_treated), " | Max treated DMS: ", np.max(all_dms_treated))
        print("Min untreated DMS: ", np.min(all_dms_untreated), " | Max untreated DMS: ", np.min(all_dms_untreated))
        print("Ratio of median: ", np.median(all_dms_treated)/np.median(all_dms_untreated))

        print('----------------')


def combine_replicates(seismic_df, min_correlation:float=0.8):

    '''

    Combine replicates A and B into one, and compute correlation between them
    Also subtract untreated signal from treated ones (bias correction)
    Only keep references that have both replicates, and good correlation

    Parameters:
        - seismic_df (pd.DataFrame): Dataframe with all the information from seismic output
        - min_cov (int): Minimum number of bases covered for a base to be used
        - min_correlation (float): Minimum correlation between A and B for a reference to be kept

    Returns:
        - combined_df (pd.DataFrame): Combined dataframe with combined replicates

    Example:
    # ref1 A is kept, and untreated is removed
    # ref2 is removed because correlation is too low
    # ref3 is removed because it only has untreated
    # ref4 is combined
    >>> seismic_df = pd.DataFrame({'sample': ['1_U', '1_A', '2_A', '2_B', '3_U', '2_A', '2_B' ], \
                                'reference': ['ref1', 'ref1', 'ref2', 'ref2', 'ref3', 'ref4', 'ref4'], \
                                'plate': [1, 1, 2, 2, 2, 2, 2], \
                                'replicate': ['Untreated', 'A', 'A', 'B', 'Untreated', 'A', 'B'], \
                                'sequence': ['AAG', 'AAG', 'GCC', 'GCC', 'ACC', 'AAGCACC', 'AAGCACC'], \
                                'sub_rate': [[UKN, 0.05 , UKN], [UKN, 0.1, UKN], [0.0 , UKN, 0.5], [0.5 , UKN, 0.0], [0.1, UKN, 0.1],     [UKN, 0.1 , UKN, UKN, 0.05, 0.25, 0.05], [UKN, 0.12 , UKN, UKN, UKN, 0.25, 0.05]], \
                                'coverage': [[UKN, 2000, UKN], [UKN, 2000, UKN], [2000, UKN, 2000], [2000, UKN, 2000], [2000, UKN, 2000], [UKN, 3000, UKN, UKN, 2000, 2000, 2000], [UKN, 1000, UKN, UKN, UKN, 2000, 2000]], \
                                'n_reads': [2000, 2000, 2000, 2000, 2000, 2000, 2000] })

    >>> combined_df = combine_replicates(seismic_df, min_correlation=0.8)

    >>> true_df = pd.DataFrame({'sample': ['1_A', ['2_A', '2_B']], \
                                'reference': ['ref1', 'ref4'], \
                                'plate': [1, 2], \
                                'sequence': ['AAG', 'AAGCACC'], \
                                'sub_rate': [[UKN, 0.1, UKN], [UKN, 0.105, UKN, UKN, 0.05, 0.25, 0.05]], \
                                'coverage': [[UKN, 2000, UKN], [UKN, 4000, UKN, UKN, 2000, 4000, 4000]], \
                                'correlation': [np.nan, 0.9940074366760487], \
                                'r2': [np.nan, 0.9815384615384616] })

    >>> pd.testing.assert_frame_equal(combined_df, true_df)

    '''
    combined_df = seismic_df.copy()
    combined_df['correlation'] = np.nan
    combined_df['r2'] = np.nan

    for ref in combined_df['reference'].unique():
        ref_df = combined_df[combined_df['reference'] == ref]
        
        # Nothing to do if only untreated replicate
        if len(ref_df[ref_df['replicate'] != 'Untreated'])==0:
            continue
        
        # If we have both replicates, we can combine them  
        if 'A' in ref_df['replicate'].values and 'B' in ref_df['replicate'].values:

            # Get signals and coverage
            signal_A = np.array(ref_df[ref_df['replicate'] == 'A']['sub_rate'].values[0])
            signal_B = np.array(ref_df[ref_df['replicate'] == 'B']['sub_rate'].values[0])
            assert len(signal_A) == len(signal_B), f"Signal not matching between replicates, {ref}"
            
            cov_A = np.array(ref_df[ref_df['replicate'] == 'A']['coverage'].values[0])
            cov_B = np.array(ref_df[ref_df['replicate'] == 'B']['coverage'].values[0])

            # Compute correlation
            mask_high_cov = (signal_A!=UKN) & (signal_B!=UKN)
            res = stats.pearsonr(signal_A[mask_high_cov], signal_B[mask_high_cov])
            r2 = r2_score(signal_A[mask_high_cov], signal_B[mask_high_cov])

            # Filter out if correlation is too low
            if res[0] < min_correlation:
                continue
                
            # Combine signals
            combined_signal = np.ones_like(signal_A)*(UKN)
            combined_cov = np.ones_like(cov_A)*(UKN)

            # When both bases have good coverage, take the average
            combined_cov[mask_high_cov] = cov_A[mask_high_cov]+cov_B[mask_high_cov]
            combined_signal[mask_high_cov] = ( signal_A[mask_high_cov]*cov_A[mask_high_cov] + signal_B[mask_high_cov]*cov_B[mask_high_cov]
                                              ) /combined_cov[mask_high_cov]

            # When only one replicate has good coverage, take its signal and not the other
            combined_signal[(signal_A!=UKN) & (signal_B==UKN)] = signal_A[(signal_A!=UKN) & (signal_B==UKN)]
            combined_signal[(signal_A==UKN) & (signal_B!=UKN)] = signal_B[(signal_A==UKN) & (signal_B!=UKN)]
            
            combined_cov[(signal_A!=UKN) & (signal_B==UKN)] = cov_A[(signal_A!=UKN) & (signal_B==UKN)]
            combined_cov[(signal_A==UKN) & (signal_B!=UKN)] = cov_B[(signal_A==UKN) & (signal_B!=UKN)]

            # Add to combined_df
            combined_df.loc[len(combined_df)] = [   ref_df['sample'].tolist(), ref_df['reference'].values[0],
                                                    ref_df['plate'].values[0], 'combined',
                                                    ref_df['sequence'].values[0], 
                                                    combined_signal.tolist(),
                                                    combined_cov.tolist(), # We don't need coverage anymore
                                                    None, # We don't need n_reads anymore
                                                    res.statistic,
                                                    r2
                                                ]
            
        # We have just one replicate, we can just copy it 
        elif 'A' in ref_df['replicate'].values or 'B' in ref_df['replicate'].values:
            combined_df.loc[len(combined_df)] = ref_df[ref_df['replicate']!='Untreated'].values[0]
            combined_df.at[combined_df.index[-1], 'replicate'] = 'combined'

    return combined_df[combined_df['replicate']=='combined'].drop(columns=['n_reads', 'replicate']).reset_index(drop=True)


# Could use better normalization, maybe with a sigmoid ?
def normalize_dms(filtered_df, percentile:float=95):

    '''

    Normalize DMS signal by the median of the top percentile

    Parameters:
        - combined_df (pd.DataFrame): Combined dataframe with combined replicates
        - percentile (float): Percentile of the top DMS signal to use for normalization

    Returns:
        - combined_df (pd.DataFrame): Combined dataframe with combined replicates, with normalized DMS signal

    Example:

    >>> combined_df = pd.DataFrame({'sample': ['s1', 's1'], \
                                'reference': ['ref1', 'ref2'], \
                                'plate': [1, 2], \
                                'sequence': ['AAGCAC', 'GAUCA'], \
                                'sub_rate': [[0.1, 0.05, UKN, UKN, 0.0, UKN], [UKN, 0.2, 0.05, UKN, 0.0]], \
                                'correlation': [None, None] })

    >>> normalized_df = normalize_dms(combined_df, percentile=95)

    >>> true_df = pd.DataFrame({'sample': ['s1', 's1'], \
                                'reference': ['ref1', 'ref2'], \
                                'plate': [1, 2], \
                                'sequence': ['AAGCAC', 'GAUCA'], \
                                'sub_rate': [[0.5, 0.25, UKN, UKN, 0.0, UKN], [UKN, 1, 0.25, UKN, 0.0]], \
                                'correlation': [None, None] })

    >>> pd.testing.assert_frame_equal(normalized_df, true_df)

    '''

    normalized_df = filtered_df.copy()

    for i, row in normalized_df.iterrows():
        
        new_dms = np.array(row['sub_rate'])
        new_dms_AC = new_dms[new_dms!=UKN]
        max_dms = np.median(new_dms_AC[new_dms_AC>=np.percentile(new_dms_AC, percentile)])
        
        if max_dms == 0: max_dms = np.max(new_dms_AC)
        if max_dms == 0: max_dms = 1

        new_dms[new_dms!=UKN] = new_dms[new_dms!=UKN]/max_dms
        new_dms[new_dms>1] = 1
        normalized_df.at[i, 'sub_rate'] = new_dms.tolist()

    return normalized_df



def plot_signal_histogram(filtered_df):

    plates = filtered_df['plate'].unique().tolist()

    fig = make_subplots(rows=len(plates), cols=1, subplot_titles=[f'Plate {plate}' for plate in plates])

    for i, plate in enumerate(plates):
        plate_df = filtered_df[filtered_df['plate']==plate]
        if not len(plate_df[plate_df['replicate']!='Untreated']):
            continue
        treated_signal = np.concatenate(plate_df[plate_df['replicate']!='Untreated']['sub_rate'].values)
        treated_signal = treated_signal[treated_signal!=UKN]
        fig.add_trace(go.Histogram(x=treated_signal[treated_signal<np.mean(treated_signal)+5*np.std(treated_signal)], 
                                   histnorm='probability', marker_color='blue', opacity=0.6, name='Treated'), row=i+1, col=1)
        if 'Untreated' in plate_df['replicate'].tolist():
            untreated_signal = np.concatenate(plate_df[plate_df['replicate']=='Untreated']['sub_rate'].values)
            untreated_signal = untreated_signal[untreated_signal!=UKN]
            fig.add_trace(go.Histogram(x=untreated_signal[untreated_signal<np.mean(untreated_signal)+5*np.std(untreated_signal)], 
                                    histnorm='probability', marker_color='red', opacity=0.6, name='Untreated'), row=i+1, col=1)
        if i!=0:
            fig.update_traces(showlegend=False, col=1, row=i+1)

    fig.update_layout(height=300*len(plates))
    fig.show()


def plot_correlations_per_plate(combined_df, mode='r2'):
    assert mode in ['correlation', 'r2'], 'mode should be either correlation or r2'
    
    plates = sorted(combined_df['plate'].unique().tolist())

    fig = make_subplots(rows=len(plates), cols=1, subplot_titles=[f'Plate {plate} | mean {mode} {np.mean(combined_df[combined_df["plate"]==plate][mode]):.2f}' for plate in plates])

    for i, plate in enumerate(plates):
        plate_df = combined_df[combined_df['plate']==plate]
        metric = plate_df[mode].to_numpy()

        fig.add_trace(go.Histogram(x=metric[~np.isnan(metric)], histnorm='probability'), row=i+1, col=1)


    fig.update_layout(height=300*len(plates), showlegend=False)
    fig.show()


def plot_histogram(df):
    all_dms = np.concatenate(df['sub_rate'].tolist())
    all_dms = all_dms[all_dms != UKN]

    fig = go.Figure()
    fig.add_trace(go.Histogram(x=all_dms, nbinsx=100))
    fig.update_layout(title_text='Distribution of DMS values', xaxis_title_text='DMS', yaxis_title_text='Count')
    fig.show()


def print_summary(df):
    for plate in np.sort(df['plate'].unique()):
        print(plate, ': ',len(df[df['plate']==plate]))

    print('Total: ', len(df))


def convert_to_rnadata(normalized_df, output_path):

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # DMS dataset format
    df_dms = normalized_df[['reference', 'sequence', 'sub_rate']].copy()
    df_dms.rename(columns={'sub_rate': 'dms_signal'}, inplace=True)
    df_dms.set_index(keys='reference', drop=True, inplace=True)
    df_dms = df_dms.T
    df_dms.to_json(os.path.join(output_path, 'dms_signal.json'), indent=2)


def plot_dms_signal(normalized_df, output_path):

    output_path = os.path.join(output_path, 'plots')

    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    for _, row in normalized_df.iterrows():

        cov = np.array(row['coverage'])
        cov[cov==UKN] = None
        dms = np.array(row['sub_rate'])
        dms[dms==UKN] = None

        seq = np.array([*row['sequence']])

        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Bar(x=np.arange(len(dms))[seq=='A'], y=dms[seq=='A'], marker_color='#636EFA', name='A signal'), secondary_y=False)
        fig.add_trace(go.Bar(x=np.arange(len(dms))[seq=='C'], y=dms[seq=='C'], marker_color='#EF553B', name='C signal'), secondary_y=False)
        fig.add_trace(go.Scatter(x=np.arange(len(dms)), y=cov, fillcolor='rgba(0,0,0,0.15)', mode='none', fill='tozeroy', name='Coverage'), secondary_y=True)

        fig.update_xaxes(title_text="Position")
        fig.update_yaxes(title_text="DMS signal")
        fig.update_yaxes(title_text="Coverage", secondary_y=True)

        plate = '' if not 'plate' in row else row['plate']+1

        fig.write_image(os.path.join(output_path, f'Plate_{plate}-ref_{row["reference"]}.png'), width=1000, height=500, scale=2)
    
    
def dump_df_to_dreem_format(seismic_df, name, path='.', metadata=None):

    if not os.path.exists(path):
        os.makedirs(path)

    if metadata is None:
        metadata = {}

    out= {**{'#Sample': name, '#Date' : datetime.today().strftime('%Y-%m-%d')}, 
            **metadata}
    for plate, df in seismic_df.groupby('plate'):
        plate = str(plate)
        for ref, df2 in df.groupby('reference'):   
            assert len(df2) == 1, 'More than one row for sample {} and reference {}'.format(plate, ref)
            df2 = df2.iloc[0]
            out[ref] = {'#samples': df2['sample'], '#plate': plate}
            out[ref]['full'] = {'sequence': df2['sequence']}
            out[ref]['full']['pop_avg'] = {
                'cov': df2['coverage'],
                'sub_rate': df2['sub_rate'],
            }
    dump_json(out, os.path.join(path, name + '.json'))
