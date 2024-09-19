#!/usr/bin/env python
# coding: utf-8

import os
import sys
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import statsmodels.stats.multitest as multi #to perform benjamini hochberg correction
import seaborn as sns
import collections 



print('######################################################################################')
print('#### Fishers exact test for classical SmileSeq sequences ©Antoni Gralak_21.08.2024####')
print('######################################################################################')

print('Setting environment...')

# To handle files, multiprocessing is available
num_cores = 1


# Set paths to data and output
data_path_input = '../output/00_readin_data/input/'
data_path_eluted = '../output/00_readin_data/'

save_path = '../output/01_fishers_exact_test/'



######################################################

try:
    os.makedirs(save_path)
except FileExistsError:
    pass

try:
    os.makedirs(save_path + 'eluted_ratios/')
except FileExistsError:
    pass

try:
    os.makedirs(save_path + 'eluted_ratios_significant/')
except FileExistsError:
    pass

#######################################################    
def process_files(files):
    translate_BCs = {'BC01': 'BC1', 'BC02': 'BC2', 'BC03': 'BC3', 'BC04': 'BC4',
                      'BC05': 'BC5', 'BC06': 'BC6_CGGATTCGAA', 'BC07': 'BC7', 'BC08': 'BC8',
                      'BC09': 'BC9', 'BC10': 'BC10', 'BC11': 'BC11', 'BC12': 'BC12'
                    }
    data_path_e = data_path_eluted
    data_path_i = data_path_input
    
    
    for i, file in enumerate(files):
        TF_name = '_'.join(file.split('_')[0:4])
        BC = file.split('_')[-1]
        BC = BC.split('.')[0]
        experiment_id = file.split('_')[-2]
        kmer = file.split('_')[-3]
        
        print(f'Processing {TF_name}, {i+1} out of {len(files)}. \n')
        
        eluted_df = pd.read_csv(os.path.join(data_path_e, file), index_col=0)
        
        # Using the SCH files as they had the deepest sequencing and thus provide maximum coverage
        input_file = f'SCH_{kmer}_{translate_BCs[BC]}.csv'
        input_df = pd.read_csv(os.path.join(data_path_i, input_file), index_col=0)
        
        # make them tidy
        input_df['status'] = 'input'
        eluted_df['status'] = 'eluted'
        
        combine = pd.concat([input_df, eluted_df])
        combine = combine.reset_index(drop=True)
        
        # Operations to make the final ratios df
        seqs = []
        for seq, df in combine.groupby('kmer'):
            # making sure that all kmers are captured
            if df.shape[0] == 2:
                seqs.append(seq)

        pvals = [None] * len(seqs)
        padj = [None] * len(seqs)
        ratio = [None] * len(seqs)

        list_of_tuples = list(zip(seqs, pvals, padj, ratio))

        result_df = pd.DataFrame(list_of_tuples, columns = ['kmer', 'pval', 'p_adjust', 'ratio_e_div_i'])

        # filter combine df for kmers which have reads in input and eluted libs

        combine = combine[combine['kmer'].isin(seqs)]
        
        #########################################################################
        # Add progress bar
        # Calculate the total number of iterations
        total_iterations = len(combine.groupby('kmer'))

        # Determine how many iterations represent 10%
        chunk_size = total_iterations // 10

        # Initialize counter and progress tracker
        counter = 0
        progress = 0

        # Print the initial empty progress bar
        print("Progress: [" + " " * 10 + "]", end="\r")
        #########################################################################

        # Perform fishers exact test and correct with benjamini hochberg. One sided fishers exact to test overrepresented kmers
        for seq, df in combine.groupby('kmer'):
            seq_e = df[df['status'] == 'eluted']['count'].values[0]
            seq_i = df[df['status'] == 'input']['count'].values[0]

            r = seq_e/seq_i

            eluted_samplesize = combine[combine['status'] == 'eluted']['count'].sum(axis=0)
            input_samplesize = combine[combine['status'] == 'input']['count'].sum(axis=0)

            rest_e = eluted_samplesize - seq_e
            rest_i = input_samplesize - seq_i

            contingency_table = np.array([[seq_e, seq_i], [rest_e, rest_i]])

            odds, p = scipy.stats.fisher_exact(contingency_table, alternative = 'greater')

            result_df.loc[result_df[result_df['kmer'] == seq].index[0]]['pval'] = p
            result_df.loc[result_df[result_df['kmer'] == seq].index[0]]['ratio_e_div_i'] =  r

            # Increment counter
            counter += 1

            # Update the progress bar every time we complete another 10% of the work
            if counter % chunk_size == 0:
                progress += 1
                # Print the updated progress bar
                sys.stdout.write(f"\rProgress: [{'=' * progress}{' ' * (10 - progress)}] {progress * 10}%")
                sys.stdout.flush()

        print('\nDone!')


        result_df["p_adjust"] = multi.multipletests(result_df["pval"], method="fdr_bh")[1]

        result_df['pval'] = pd.to_numeric(result_df['pval'])
        result_df['p_adjust'] = pd.to_numeric(result_df['p_adjust'])
        result_df['ratio_e_div_i'] = pd.to_numeric(result_df['ratio_e_div_i'])
        
        
        
        # Now see how many kmers are significant
        
        threshold = 0.05

        result_df['below_threshold'] = result_df['p_adjust'] < threshold
        
        result_df.to_csv(save_path + f'eluted_ratios/{TF_name}_{kmer}_{experiment_id}_{BC}.csv')
        
        
        significant = result_df[result_df['below_threshold'] == True]
        
        if len(significant) >= 5:
            significant.to_csv(save_path + f'eluted_ratios_significant/{TF_name}_{kmer}_{experiment_id}_{BC}.csv')
            #verdict = 'ENRICHMENT'
        #else:
            #verdict = 'FAIL'
        
        
        #if verdict == 'ENRICHMENT':
        #    significant.to_csv(f'eluted_ratios_significant/{TF_name}_{kmer}_{experiment_id}_{BC}.csv')
        #    sns.swarmplot(data=significant.nsmallest(100, columns='p_adjust'), x="ratio_e_div_i", rasterized=True)
        #    plt.savefig(f'eluted_ratios_swarmplot/{TF_name}_{kmer}_{BC}_swarm.pdf')
        #    plt.close()
            
    # Get summary stats for each processed TF and plot if needed
    
    #summary = pd.concat([summary] + summary_dfs_list, ignore_index=True)
    
    #return summary

#######################################################  


# Run script with n cores to extract kmers of all eluted reads

if __name__ == "__main__":
    
    files = os.listdir(data_path_eluted)
    
    num_processes = num_cores
    
    # If needed to select a certain kmer.
    f_files = [file for file in files if '9mer' in file]
    
    #chunks = [files[i::num_processes] for i in range(num_processes)]
    chunks = [f_files[i::num_processes] for i in range(num_processes)]
    
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        
        pool.map(process_files, chunks)
       