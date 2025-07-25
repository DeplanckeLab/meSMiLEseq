#!/usr/bin/env python
# coding: utf-8
# %%

# %%


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
print('#### Fishers exact test with SmileSeq sequences ©Antoni Gralak_20.08.24##')
print('#########################################################################')
print('Setting env...')
this_path = os.getcwd()
sys.path.append(this_path)

import functions

# Load metadata
metadata = pd.read_csv('../exemplary_data/metadata.csv', sep=';')

data_path = '../output/01_kmer_analysis/'

save_path = '../output/02_fishers_exact_test/'

try:
    os.mkdir(save_path)
except FileExistsError:
    pass



# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in sequences generated in 01_kmer_analysis and calculates pvalues
based on fishers exact test (assuming a hypergeometric distribution). Checks for enriched kmers (calculates the probability
of obtaining the same amount or more same kmers). Pvalues are corrected via benjamini-hochberg).""")

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


# iterate over the TFs and kmers in SmS and load the file

print('Starting p_val calculations.')
for TF in to_be_analyzed:
    print(f'Starting with {TF}.')
    for kmer in kmers:
        
        large_df = pd.read_csv(data_path + f'{experiment_name}_{TF}_{kmer}mer_enrichment.csv')

        seqs = []
        mods = []
        for seq, df in large_df.groupby('kmer'):
            if df[df['status'] == 'input'].shape[0] == df[df['status'] == 'eluted'].shape[0]:
                if (df[df['mod'] == 'nonmethl'].shape[0] == 2) | (df[df['mod'] == 'methl'].shape[0] == 2):
                    modifications = list(df[df['status'] == 'eluted']['mod'].values)
                    seq = [seq] * len(modifications)


                    seqs.extend(seq)
                    mods.extend(modifications)
        
        pvals = [None] * len(seqs)
        padj = [None] * len(seqs)

        list_of_tuples = list(zip(seqs, mods, pvals, padj))

        result_df = pd.DataFrame(list_of_tuples, columns = ['kmer', 'mod', 'pval', 'p_adjust'])

        # Filter dataframe for kmers in result_df, iterate over each kmer, create contingency table and perform fisher exact test


        large_df = large_df[large_df['kmer'].isin(seqs)]


        for seq, df in large_df.groupby('kmer'):
            for mod in ['methl', 'nonmethl']:
                if df[df['mod'] == mod].empty:
                    continue
                else:
                    seq_e = df[(df['status'] == 'eluted') & (df['mod'] == mod)]['count'].values[0]
                    seq_i = df[(df['status'] == 'input') & (df['mod'] == mod)]['count'].values[0]
        
                    eluted_samplesize = large_df[(large_df['status'] == 'eluted') & (large_df['mod'] == mod)]['count'].sum(axis=0)
                    input_samplesize = large_df[(large_df['status'] == 'input') & (large_df['mod'] == mod)]['count'].sum(axis=0)
        
                    rest_e = eluted_samplesize - seq_e
                    rest_i = input_samplesize - seq_i
        
                    contingency_table = np.array([[seq_e, seq_i], [rest_e, rest_i]])

                    odds, p = scipy.stats.fisher_exact(contingency_table, alternative = 'greater')
        
                    result_df.loc[result_df[(result_df['kmer'] == seq) & (result_df['mod'] == mod)].index[0]]['pval'] = p
        

        result_df["p_adjust"] = multi.multipletests(result_df["pval"], method="fdr_bh")[1]



        print('Done!')

        result_df.to_csv(save_path + f'{experiment_name}_{TF}_{kmer}mer_pvalues.csv', index=False)
