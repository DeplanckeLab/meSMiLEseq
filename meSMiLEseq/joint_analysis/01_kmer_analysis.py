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
sys.path.append(this_path)


# %%


import my_functions


# %%


# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in sequences generated in 00_read_in and splits the random region
of 24 nucleotides into kmers k=[6, 7, 8, 9]. It also calculates a p_value for all kmers 
(multiple test correction benjamini-hochberg).""")

parser.add_argument('-sms', '--sms_name', type=str, help='Smile-seq experiment number. E.g. SmSAG01.', required=True)
parser.add_argument('-k', '--kmer', type=int, nargs='+', help="""kmer size. By default k=[6, 7, 8, 9]. Parse multiple
 or a single integer.""")
parser.add_argument('-bc', '--Barcode', type=str, nargs='+', help="""Barcodes to be included. By default all BC1 to 12.
Parse multiple or single BCs.""")
#parser.add_argument('-gr', '--graphs', type=str, help='Do you want summary graphs? True/False', required=True)

# Parse the command line arguments
args = parser.parse_args()
arguments = vars(args)


# %%


experiment_name = arguments['sms_name'] 


# %%


if experiment_name in ['SmSAG01', 'SmSAG02', 'SmSAG03']:
    input_data_path = '/home/gralak/updepla/users/gralak/NAS2/SmileSeq_paper/SmileSeq_experiments/inputs/20220315_input/00_read_in_data/output/'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"

elif experiment_name in ['SmSAG04', 'SmSAG05', 'SmSAG06', 'SmSAG07', 'SmSAG08']:
    input_data_path = '/home/gralak/updepla/users/gralak/NAS2/SmileSeq_paper/SmileSeq_experiments/inputs/20220915_input/00_read_in_data/output/'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"

elif experiment_name in ['SmSAG09', 'SmSAG10', 'SmSAG11', 'SmSAG12', 'SmSAG13', 'SmSAG14']:
    input_data_path = '/home/gralak/updepla/users/gralak/NAS2/SmileSeq_paper/SmileSeq_experiments/inputs/20230228_input/00_read_in_data/output/'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"

elif experiment_name in ['SmSAG15', 'SmSAG16', 'SmSAG17', 'SmSAG18', 'SmSAG19', 'SmSAG20', 'SmSAG21', 'SmSAG22', 'SmSAG23']:
    input_data_path = '/home/gralak/updepla/users/gralak/NAS2/SmileSeq_paper/SmileSeq_experiments/inputs/20230503_input/00_read_in_data/output/'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"

else:
    print("-sms needs to be a SmileSeq name such as SmSAG01 (possible options 01 to 23). Stopping script.")
    sys.exit(1)


# %%


data_path = '/home/gralak/updepla/users/gralak/NAS2/SmileSeq_paper/SmileSeq_experiments/' + experiment_name + '/00_read_in_data/output/'

save_path = '/home/gralak/updepla/users/gralak/SmileSeq_paper/meSMiLEseq_joint_analysis/' + experiment_name + '/'


# %%


try:
    os.makedirs(save_path + '01_kmer_analysis/')
except FileExistsError:
    pass

try:
    os.makedirs(save_path + '02_fishers_exact_test/')
except FileExistsError:
    pass

output_path_01 = save_path + '01_kmer_analysis/'
output_path_02 = save_path + '02_fishers_exact_test/'



# %%


if experiment_name == 'SmSAG01':
    barcodes = pd.read_csv(this_path + '/' + experiment_name + '/Barcodes.csv', sep='\t', header=None)
    TFs = pd.read_csv(this_path + '/' + experiment_name + '/analysed_TFs.csv', sep=';')
    
    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"
    
elif experiment_name in ['SmSAG02', 'SmSAG03']:
    barcodes = pd.read_csv(this_path + '/' + experiment_name + '/Barcodes.csv', sep=',', header=None)
    TFs = pd.read_csv(this_path + '/' + experiment_name + '/analysed_TFs.csv', sep=';')

    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"
    

else:
    barcodes = pd.read_csv(this_path + '/' + experiment_name + '/Barcodes.csv', header=None)
    TFs = pd.read_csv(this_path + '/' + experiment_name + '/analysed_TFs.csv')
    
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"


TFs = TFs[barcodes.loc[:,2]]

#reindex
TFs['index'] = np.arange(0, len(TFs), step=1)
TFs = TFs.rename(index=TFs['index'])


# %%


names = barcodes[barcodes[2] == True][0].values.tolist()


# %%


# Define what will be analyzed

if arguments['Barcode']:
    to_be_analyzed = arguments['Barcode']
    
else:
    to_be_analyzed = ['BC1','BC2','BC3','BC4','BC5','BC6','BC7','BC8','BC9','BC10','BC11','BC12']


if arguments['kmer']:
    kmers = arguments['kmer']
else:
    kmers = [6, 7, 8, 9]


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


print(f'Starting analysis for {experiment_name}.')
for iterator, BC in enumerate(to_be_analyzed):
    
    print('loading necessary datasets...')
    print(f'loading {experiment_name}, {BC}...')
    input_df = pd.read_csv(input_data_path + f'{BC}_contamination_filtered.csv')
    eluted_df = pd.read_csv(data_path + f'{BC}_contamination_filtered.csv')



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


    name = TFs[barcodes[0] == BC]['Proteins'].values[0]

############################################################################################
# split into kmers

    for kmer in kmers:
        print(f'Starting with {kmer}mer for {name}...')
        large_df = my_functions.kmer_counting(df=[input_df, eluted_df],
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


        large_df.to_csv(output_path_01 + f'{name}_{str(kmer)}mer_enrichment.csv', index=False)
        print(str(kmer) + f'mer dataframe for {name} saved!')



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

        result_df.to_csv(output_path_02 + f'{name}_{str(kmer)}kmer_pvalues.csv', index=False)
        print('Done!')
        print('#' * 10)


        
print('Finished!')

