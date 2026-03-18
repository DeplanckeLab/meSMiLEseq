#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import argparse

import os
import sys

import numpy as np
import pandas as pd

import scipy
import statsmodels.stats.multitest as multi #to perform benjamini hochberg correction




# %%


print('#########################################################################')
print('#### Data Wrangling with SmileSeq sequences ©Antoni Gralak_19.05.2025####')
print('#########################################################################')
print('Setting env...')
this_path = os.getcwd()
sys.path.append(os.path.join(this_path, '..'))
import utils




# Load metadata
metadata = pd.read_csv('../metadata.csv')

data_path = '../output_joint_analysis/00_read_in_data/'

save_path = '../output_joint_analysis/'
# %%

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
#parser.add_argument('-gr', '--graphs', type=str, help='Do you want summary graphs? True/False', required=True)

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




output_path_01 = save_path + '01_kmer_analysis/'
output_path_02 = save_path + '02_fishers_exact_test/'

os.makedirs(output_path_01, exist_ok=True)
os.makedirs(output_path_02, exist_ok=True)




# %%


# Since I some barcodes might have mutations, these two functions will allow to identify the mBCs with
# a Hamming distance 2 and chage them to the respective BC.


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
        if hamming_distance(input_str, s) < 2:
            similar_strings.append(True)
        else:
            similar_strings.append(False)
    return similar_strings


# ### For loop over different kmer lengths -> save everything


# %%


print('Starting analysis.')
for TF in to_be_analyzed:
    
    print(f'loading necessary datasets... {TF}')

    position_on_chip = metadata[(metadata['experiment'] == experiment_name) & (metadata['TF'] == TF)]['Chip_pos'].values[0]
    
    input_df = pd.read_csv(data_path + f'{input_id}_{position_on_chip}_raw_data.csv')
    eluted_df = pd.read_csv(data_path + f'{experiment_name}_{TF}_raw_data.csv')



    ##########################################################################################
    # input first, change all mBCs with hamming dsitance 1 or less into the corresponding mBCs.
    mBCs = input_df['methl']
    
    #change all mBC which are hamming distance 1 from AGTA to AGTA
    similar_strings = find_similar_strings(methylated_BC, mBCs)
    
    input_df.loc[similar_strings, 'methl'] = methylated_BC
    
    #change all mBC which are hamming distance 1 from GAAT to GAAT
    similar_strings = find_similar_strings(unmethylated_BC, mBCs)
    
    input_df.loc[similar_strings, 'methl'] = unmethylated_BC
    
    
    
    # the same for eluted
    mBCs = eluted_df['methl']
    
    #change all mBC which are hamming distance 1 from AGTA to AGTA
    similar_strings = find_similar_strings(methylated_BC, mBCs)
    
    eluted_df.loc[similar_strings, 'methl'] = methylated_BC
    
    #change all mBC which are hamming distance 1 from GAAT to GAAT
    similar_strings = find_similar_strings(unmethylated_BC, mBCs)
    
    eluted_df.loc[similar_strings, 'methl'] = unmethylated_BC


    print('Data loaded.')

############################################################################################
# split into kmers

    for kmer in kmers:
        print(f'Starting with {kmer}mer for {TF}...')
        large_df = utils.kmer_counting(df=[input_df, eluted_df],
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


        large_df.to_csv(output_path_01 + f'{experiment_name}_{TF}_{kmer}mer_enrichment.csv', index=False)
        print(str(kmer) + f'mer dataframe for {TF} saved!')



        ############################################################################################
        # Continuation with pvalues
        print('#' * 10)
        print('Continuing with p-value calculation...')
        print("Disclaimer: this part can be found separately in '02_fishers_exact_test.py'. For convenience it is included in this script.")
        
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
        
                    result_df.loc[result_df[(result_df['kmer'] == seq) & (result_df['mod'] == mod)].index[0], 'pval'] = p
        

        result_df["p_adjust"] = multi.multipletests(result_df["pval"], method="fdr_bh")[1]



        print('Done!')

        result_df.to_csv(output_path_02 + f'{experiment_name}_{TF}_{kmer}mer_pvalues.csv', index=False)
        print('Done!')
        print('#' * 10)


        
print('Finished!')

