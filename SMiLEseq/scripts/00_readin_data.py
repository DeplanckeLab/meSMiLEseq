#!/usr/bin/env python
# coding: utf-8


print('#################################################################################')
print('#### Reading data from classical SmileSeq sequences ©Antoni Gralak_21.08.2024####')
print('#################################################################################')

print('Setting environment...')

import os
import sys
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import collections

sys.path.append('../../meSMiLEseq/scripts/')

import functions

# To handle eluted files, multiprocessing is available
num_cores = 1

# First define parameters to process input files

BCs = ['BC7'] # adjust for BC7 (this script is merely an example)
kmers = [9] # adjust for kmers if needed

# Set paths to data and output
data_path_input = '../data/input/'
data_path_eluted = '../data/eluted_raw/'
save_path_input = '../output/00_readin_data/input/'
save_path = '../output/00_readin_data/'

######################################################
try:
    os.makedirs(save_path_input)
except FileExistsError:
    pass


try:
    os.makedirs(save_path)
except FileExistsError:
    pass

# These functions are customized for classical SMiLEseq data

def kmer_counting(df, kmer=6):
    """Here, the random region from the classical SMS (40 bp) are divided into kmers (default 6) and counted."""
    ks = [df[0][i][x:(x+kmer)] for i in range(len(df)) for x in range(40-kmer+1)] #generates all kmers
    
    ks = pd.DataFrame.from_dict(dict(collections.Counter(ks)), orient='Index')
    ks = ks.reset_index()
    
    ks = ks.rename(columns={0:'count', 'index':'kmer'})
    return ks


def ecdf(data):
    x = np.sort(data)
    y = np.arange(1, len(data) + 1) / len(data)
    return x, y

print('Starting with inputs...')
for BC in BCs:
    
    two = f'SCH_INP_Trimmed_{BC}.fastq.gz'
    
    two2 = functions.readin_fastq(data_path_input,two)
    two2 = two2.reset_index(drop=True)

    
    for kmer in kmers:

        two22 = kmer_counting(two2, kmer=kmer)
        two22.to_csv(save_path_input + f'SCH_{kmer}mer_{BC}.csv')


        #######################################
        # Z-transform

        mu2 = np.mean(two22['count'])
        std2 = np.std(two22['count'])
        two22['zscore'] = two22['count'].apply(lambda x : (x-mu2)/std2)

        #######################################

        x2 = np.array(two22['zscore'])
        x2, y2 = ecdf(x2)

        #######################################
        #plotting
        fig, ax = plt.subplots(figsize=(6,6))


        #input data 2

        ax.plot(x2, y2, marker=None, linestyle='-', color='black', linewidth=2)

        ax.grid(True, alpha=0.2)
        plt.savefig(save_path_input + f'{BC}_{kmer}mer_ecdf_inputs.pdf')
        plt.show()
        print(f'Done with {BC}!')

print('Finished! Now eluted...')

# Next, process eluted files
def process_files(files):
    """Pass a list of files that need to be processed."""
    raw_data_path = data_path_eluted
    kmers = [9]#[6,7,8]
    plot = False
    
    for i, file in enumerate(files):
        TF_name = '_'.join(file.split('_')[0:4])
        BC = file.split('_')[-1]
        BC = BC.split('.')[0]
        experiment_id = file.split('_')[-2]
        print(f'Working on {TF_name}, TF {i+1} out of {len(files)}. \n')
        
        df = functions.readin_fastq(raw_data_path,file)
        df = df.reset_index(drop=True)
        
        for kmer in kmers:
            # calc kmers
            df_kmer = kmer_counting(df, kmer=kmer)
            df_kmer.to_csv(save_path + f'{TF_name}_{kmer}mer_{experiment_id}_{BC}.csv')
            
            if plot:
                # ecdf
                mu1 = np.mean(df_kmer['count'])
                std1 = np.std(df_kmer['count'])
                df_kmer['zscore'] = df_kmer['count'].apply(lambda x : (x-mu1)/std1)

                x1 = np.array(df_kmer['zscore'])
                x1, y1 = ecdf(x1)

                fig, ax = plt.subplots(figsize=(6,6))

                #input data 1 

                ax.plot(x1, y1, marker=None, linestyle='-', color='green', linewidth=2)
                ax.grid(True, alpha=0.2)
                plt.savefig(save_path + f'{TF_name}_{kmer}mer_{experiment_id}_{BC}_ecdf.pdf')
                plt.close()

# Run script with n cores to extract kmers of all eluted reads

if __name__ == "__main__":
    
    files = os.listdir(data_path_eluted)
    
    num_processes = num_cores
    
    chunks = [files[i::num_processes] for i in range(num_processes)]
    
    # Create a pool of processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Map the process_files function to each chunk
        pool.map(process_files, chunks)
    
    print('Done!')