#!/usr/bin/env python
# coding: utf-8

import os
import sys
import multiprocessing
import re
import numpy as np
import pandas as pd

sys.path.append('../../meSMiLEseq/scripts/')

import functions


print('#####################################################################################################################')
print('#### Extracting significant fastqs for motif discovery of classical SmileSeq sequences ©Antoni Gralak_21.08.2024 ####')
print('#####################################################################################################################')

print('Setting environment...')

# To handle files, multiprocessing is available
num_cores = 1

# Set paths to data and output
raw_eluted_p = '../data/eluted_raw/'
    
significant_kmers_p = '../output/01_fishers_exact_test/eluted_ratios_significant/'

output_path_ = '../output/02_extract_significant_fastqs/'



######################################################

try:
    os.makedirs(output_path_)
except FileExistsError:
    pass

######################################################


def choose_kmer(kmer_files):
    'Selects the longest kmer (9mer) to filter fastqs.'
    for i in kmer_files:
        if '9mer' in i:
            return i, '9mer'

    for j in kmer_files:
        if '8mer' in j:
            return j, '8mer'
    
    for k in kmer_files:
        if '7mer' in k:
            return k, '7mer'
    
    for l in kmer_files:
        if '6mer' in l:
            return l, '6mer'


# Function to check if any snippet exists in a given text


def check_snippet(text, snippets):
    for snippet in snippets:
        if snippet in text:
            return True
    return False


def process_files(significant_ids):
    # Set env
    raw_eluted_path = raw_eluted_p
    
    significant_kmers = significant_kmers_p 

    output_path = output_path_
    
    
    
    # all raw fastq files and all significant kmer files
    all_TF = os.listdir(raw_eluted_path)
    all_kmers = os.listdir(significant_kmers)
    
    
    # take significant ID, search for significant kmers and filter the raw fastq for reads with sig kmers
    for UT_ID in significant_ids:
        
        
        # create regex pattern with UT identifier and search for it in eluted fastqs and kmers.
        pattern = re.compile(rf'{re.escape(UT_ID)}.*')
        
        # Here are all the raw fastqs which should be processed, saved as list
        raw_fastqs = [ID for ID in all_TF if pattern.match(ID)]
        
        # Same process to get significant kmer files
        kmer_files = [ID for ID in all_kmers if pattern.match(ID)]
        
        # Select the largest kmer (computed up to 8mers) and readin df
        chosen_kmer, length = choose_kmer(kmer_files)
        chosen_kmer_df = pd.read_csv(os.path.join(significant_kmers, chosen_kmer), index_col=0)
        
        # readin fastqs and filter rows for those with significant kmer and save it all
        for file in raw_fastqs:
            name = file.split('.')[0]
            raw_fastq = functions.readin_fastq(raw_eluted_path,file)
            raw_fastq = raw_fastq.reset_index(drop=True)
            
            # filter out rows which do not contain significant kmers
            filtered_fastq = raw_fastq[raw_fastq[0].apply(lambda x: check_snippet(x, chosen_kmer_df['kmer']))]
            filtered_fastq.to_csv(os.path.join(output_path, f'{name}_{length}_filtered.fastq'))
            
       
# Run script and save under eluted_significant_raw

if __name__ == "__main__":
    
    files = os.listdir(significant_kmers_p)
    ids = ['_'.join(name.split('_')[0:2]) for name in files]
    # remove duplicated UT_IDs
    significant_ids = list(set(ids))
    
    num_processes = num_cores
    
    chunks = [significant_ids[i::num_processes] for i in range(num_processes)]
    
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Call process_files function for each chunk
        pool.map(process_files, chunks)