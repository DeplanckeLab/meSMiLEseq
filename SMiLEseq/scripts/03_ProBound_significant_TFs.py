#!/usr/bin/env python
# coding: utf-8

import pyProBound_operator as pbo
import os
import multiprocessing
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import re
import sys

sys.path.append('../../meSMiLEseq/scripts/')

import functions

print('###################################################################################')
print('#### Motif discovery of classical SmileSeq sequences ©Antoni Gralak_21.08.2024 ####')
print('###################################################################################')

print('Setting environment...')

# Run script with n cores to create ProBound model, careful ProBound takes 4 cores
num_cores = 1


# Set env and reused variables, paths to data and output
calc_path_p = '../output/03_ProBound_significant_TFs/calculations/'
psam_path_p = '../output/03_ProBound_significant_TFs/psams/'
    
# Change paths if needed! E.g. take subset or 9mers
raw_eluted_p = '../output/02_extract_significant_fastqs/'
raw_input_p = '../data/input/'



######################################################

try:
    os.makedirs(calc_path_p)
except FileExistsError:
    pass

try:
    os.makedirs(psam_path_p)
except FileExistsError:
    pass

######################################################


# Use filtered fastqs to infer motif.

def process_files(files):
    
    # Set env and reused variables
    calc_path = calc_path_p
    psam_path = psam_path_p
    
    # Change paths if needed! E.g. take subset or 9mers
    raw_eluted_path = raw_eluted_p
    raw_input_path = raw_input_p
    
    
    translate_BCs = {'BC01': 'BC1', 'BC02': 'BC2', 'BC03': 'BC3', 'BC04': 'BC4',
                    'BC05': 'BC5', 'BC06': 'BC6_CGGATTCGAA', 'BC07': 'BC7', 'BC08': 'BC8',
                    'BC09': 'BC9', 'BC10': 'BC10', 'BC11': 'BC11', 'BC12': 'BC12'
                    }
    
    BC_to_sequence = {'BC01': 'CGTATGAATC', 'BC02': 'AAGTTAGAAG', 'BC03': 'TTATGGAACA', 'BC04': 'GTTACACGCC',
                    'BC05': 'CGGATTAGAC', 'BC06': 'AGAAATGTGG', 'BC07': 'AACGAGTGTG', 'BC08': 'AGTTAAAGTC',
                    'BC09': 'CATCGTGGTC', 'BC10': 'CTCAGATATC', 'BC11': 'GACGATGAGA', 'BC12': 'TCAGGCCACA'
                     }
    
    binding_mode_size = [15]# [9, 15, 25]
    
    for file in files:
        
        # Get specs
        UT_ID = '_'.join(file.split('_')[0:2])
        TF_name = file.split('_')[2]
        experiment_id = file.split('_')[-4]
        BC = file.split('_')[-3]
        filtered = file.split('_')[-2]
        
        # readin eluted fastq (the fastqs aren't real fastqs, they are csvs) 
        fastq = pd.read_csv(os.path.join(raw_eluted_path,file), index_col=0)
        fastq.columns = ['sequence']
        
        
        # readin input fastq
        input_fastq = functions.readin_fastq(raw_input_path, f'SCH_INP_Trimmed_{translate_BCs[BC]}.fastq.gz')
        input_fastq = input_fastq.reset_index(drop=True)
        input_fastq.columns = ['sequence']
        
        # filter out Ns
        input_fastq = input_fastq[~input_fastq["sequence"].str.contains("N")]
        fastq = fastq[~fastq["sequence"].str.contains("N")]
        
        
        # set flanks for ProBound
        left_flank = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" + BC_to_sequence[BC]
        right_flank = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
        
        
        # create output path to write models and psams
        for binding_mode in binding_mode_size:
            
            TF_path_calc = os.path.join(calc_path, TF_name, '_'.join([UT_ID,experiment_id,filtered]), f'{binding_mode}bps')
            TF_path_psam = os.path.join(psam_path, TF_name, '_'.join([UT_ID,experiment_id,filtered]), f'{binding_mode}bps')
            try:
                os.makedirs(TF_path_calc)
            except FileExistsError:
                pass
            try:
                os.makedirs(TF_path_psam)
            except FileExistsError:
                pass
            
            #########################################
            # Create config for ProBound and set env#
            #########################################
            print('Create config for ProBound and set env')

            outputfile = os.path.join(TF_path_calc, "f_output.tsv.gz")
            
            count_table = pbo.build_count_table(input_fastq, fastq,
                                    output_filename=outputfile, gzip=True)
            
            
            
            # the default tested configuration for smile seq with three binding modes
            config = pbo.generate_SMiLE_seq_configuration(outputfile,
                                                        variable_region_length=40,
                                                        left_flank=left_flank,
                                                        right_flank=right_flank,
                                                        binding_mode_flank=30, 
                                                          # this must be smaller than the left and right flank size
                                                        binding_modes=2,
                                                        binding_mode_size=binding_mode)
            
            
            
            basename=TF_name + '_testmodel'
            # this is the identification of the new model (<base_name>.models.json)

            config.alter_output(output_path=TF_path_calc, 
                                base_name=basename, 
                                print_trajectory=True, 
                                # if true, generates several files in the output path, one set for each binding mode
                                # <base_name>.trajectory.component<binding mode index>-<desc of file>.csv
                                verbose=False) # flipping this switch does not seem to do very much tbh
            
            # Once you are done with the config modifications, write it to file. 
            # For ProBound, only what is written to the config file counts!
            
            config_filename = os.path.join(TF_path_calc, TF_name + "_config.json")
            config.print_json(config_filename)
            
            # Run ProBound
            print(f'Running ProBound for {TF_name}, binding size {binding_mode}, {filtered} filtered.')
            pbo.run_probound(config_filename, 
                             full_config_file=os.path.join(TF_path_calc, "tmp.fullconfig.json"),
                             save_output=os.path.join(TF_path_calc, "tmp.optimization.out"),
                             cleanup_verbose=True)
            
            # Retrieve psam
            result_filename = pbo.get_psam(os.path.join(TF_path_calc, f"{basename}.models.json"))
            
            # Extract psam and plot dGG matrix
            for j, psam in enumerate(result_filename):
                
                # If binding mode is all but negative, don't even plot, just skip
                if (psam < 0).all().all():
                    psam.to_csv(os.path.join(TF_path_psam, TF_name + f"_bindingmode_{str(j + 1)}_BELOW_ZERO.csv"))
                    continue
                else:
                    psam.to_csv(os.path.join(TF_path_psam, TF_name + f"_bindingmode_{str(j + 1)}.csv"))

                    fig, ax = plt.subplots(1,1,figsize=[10,6])
                    logo = logomaker.Logo(result_filename[j],
                                        shade_below=0.5,
                                        ax=ax,
                                        fade_below=0.5,
                                        color_scheme={'A':'#66a61e', 'C':'#7570b3','G':'#ffc809','T':'#d95f02','m':'#a6cee3'}
                                        )
                    # style using Logo methods
                    logo.style_spines(visible=False)
                    logo.style_spines(spines=['left', 'bottom'], visible=True)
                    logo.style_xticks(rotation=90, fmt='%d', anchor=0)

                    # style using Axes methods
                    logo.ax.set_ylabel("$-\Delta \Delta G$ (kcal/mol)", labelpad=-1)
                    logo.ax.xaxis.set_ticks_position('none')
                    logo.ax.xaxis.set_tick_params(pad=-1)
                    #logo.ax.set_ylim([-6, 4])

                    fig.suptitle(f"{TF_name} bindingmode {str(j + 1)}")

                    fig.savefig(os.path.join(TF_path_psam, TF_name + f"_bindingmode_{str(j + 1)}_logo.pdf"), format='pdf')
                    fig.savefig(os.path.join(TF_path_psam, TF_name + f"_bindingmode_{str(j + 1)}_logo.png"), format='png')
                    plt.close()
            print(f'Done with {TF_name}, binding size {binding_mode}, {filtered} filtered.')    
# Run script with n cores to create ProBound model, careful ProBound takes 4 cores

if __name__ == "__main__":
    
    all_TF = os.listdir(raw_eluted_p)
    
    num_processes = num_cores
    
    chunks = [all_TF[i::num_processes] for i in range(num_processes)]
    
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Call process_files function for each chunk
        pool.map(process_files, chunks)
