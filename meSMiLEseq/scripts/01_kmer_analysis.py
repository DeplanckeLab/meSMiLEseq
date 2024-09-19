#!/usr/bin/env python
# coding: utf-8

import argparse

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import statsmodels.stats.multitest as multi #to perform benjamini hochberg correction

from collections import Counter


print('#########################################################################')
print('#### k-mer analysis with SmileSeq sequences ©Antoni Gralak_19.08.2024####')
print('#########################################################################')
print('Setting env...')
this_path = os.getcwd()
sys.path.append(this_path)

import functions

# Load metadata
metadata = pd.read_csv('../exemplary_data/metadata.csv', sep=';')

data_path = '../output/00_read_in_data/'

save_path = '../output/01_kmer_analysis/'

try:
    os.mkdir(save_path)
except FileExistsError:
    pass

# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in sequences generated in 00_read_in and splits the random region
of 24 nucleotides into kmers k=[6, 7, 8, 9]. It also calculates a p_value for all kmers 
(multiple test correction benjamini-hochberg).""")

parser.add_argument('-sms', '--sms_name', type=str, help='Smile-seq experiment number. E.g. exp1.', required=True)
parser.add_argument('-k', '--kmer', type=int, nargs='+', help="""kmer size. By default k=[6, 7, 8, 9]. Parse multiple
 or a single integer.""")
parser.add_argument('-tf', '--Transcription_factor', type=str, nargs='+', help="""TFs to be included. By default all TFs
 that were approved in the experiment. Parse multiple or single TFs.""")


# Parse the command line arguments
args = parser.parse_args()
arguments = vars(args)


experiment_name = arguments['sms_name']


if experiment_name in ['exp1', 'exp2', 'exp3']:
    input_id = 'input1' 
    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"
elif experiment_name in ['exp4', 'exp5', 'exp6', 'exp7', 'exp8']:
    input_id = 'input2'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
elif experiment_name in ['exp9', 'exp10', 'exp11', 'exp12', 'exp13', 'exp14']:
    input_id = 'input3'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
elif experiment_name in ['exp15', 'exp16', 'exp17', 'exp18', 'exp19', 'exp20', 'exp21', 'exp22', 'exp23']:
    input_id = 'input4'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
else:
    print("-sms needs to be a experiment ID, e.g. exp1 (possible options 1 to 23). Stopping script.")
    sys.exit(1)


# Define what will be analyzed

if arguments['Transcription_factor']:
    to_be_analyzed = arguments['Transcription_factor']
    
else:
    to_be_analyzed = list(metadata[(metadata['experiment'] == experiment_name) & (metadata['approved'] == True)]['TF'])


if arguments['kmer']:
    kmers = arguments['kmer']
else:
    kmers = [6, 7, 8, 9]




# Since I some barcodes might have mutations, these two functions will allow to identify the mBCs with
# a Hamming distance 2 and change them to the respective BC.


def hamming_distance(s1, s2):
    """
    Calculates the Hamming distance between two strings s1 and s2.
    """
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_similar_strings(input_str, strings):
    """
    Finds strings in the list `strings` that have a Hamming distance of up to 2
    from the input string `input_str`.
    """
    similar_strings = []
    for s in strings:
        if hamming_distance(input_str, s) <= 2:
            similar_strings.append(True)
        else:
            similar_strings.append(False)
    return similar_strings


# Loop over different kmer lengths -> save everything





print('Starting analysis.')
for TF in to_be_analyzed:
    
    print(f'loading necessary datasets... {TF}')
    
    position_on_chip = metadata[(metadata['experiment'] == experiment_name) & (metadata['TF'] == TF)]['Chip_pos'].values[0]
    input_df = pd.read_csv(data_path + f'{input_id}_{position_on_chip}_raw_data.csv')
    eluted_df = pd.read_csv(data_path + f'{experiment_name}_{TF}_raw_data.csv')



    ##########################################################################################
    # input first, change all mBCs with hamming dsitance 2 or less into the corresponding mBCs.
    mBCs = input_df['methl']
    
    #change all mBC which are hamming distance 2 from AGTA to AGTA
    similar_strings = find_similar_strings(methylated_BC, mBCs)
    
    input_df.loc[similar_strings, 'methl'] = methylated_BC
    
    #change all mBC which are hamming distance 2 from GAAT to GAAT
    similar_strings = find_similar_strings(unmethylated_BC, mBCs)
    
    input_df.loc[similar_strings, 'methl'] = unmethylated_BC
    
    
    
    # the same for eluted
    mBCs = eluted_df['methl']
    
    #change all mBC which are hamming distance 2 from AGTA to AGTA
    similar_strings = find_similar_strings(methylated_BC, mBCs)
    
    eluted_df.loc[similar_strings, 'methl'] = methylated_BC
    
    #change all mBC which are hamming distance 2 from GAAT to GAAT
    similar_strings = find_similar_strings(unmethylated_BC, mBCs)
    
    eluted_df.loc[similar_strings, 'methl'] = unmethylated_BC


    print('Data loaded.')
############################################################################################
# split into kmers

    for kmer in kmers:
        print(f'Starting with {kmer}mer for {TF}...')
        large_df = functions.kmer_counting(df=[input_df, eluted_df],
                                              kmer=kmer,
                                              status=['input','eluted'],
                                              one_df=False,
                                              mbc=[methylated_BC, unmethylated_BC]
                                             )
        
        #Add info about CpG presence in kmer
        CpG = []
        for i in range(len(large_df)):
            if 'CG' in large_df['kmer'].values[i]:
                CpG.append('present')
            else:
                CpG.append('not present')

        large_df.insert(4,'CpG', CpG)


        large_df.to_csv(save_path + f'{experiment_name}_{TF}_{kmer}mer_enrichment.csv', index=False)
        print(str(kmer) + f'mer dataframe for {TF} saved!')



        
print('Finished!')

