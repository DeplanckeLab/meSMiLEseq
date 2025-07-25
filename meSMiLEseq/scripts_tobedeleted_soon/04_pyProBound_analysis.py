#!/usr/bin/env python
# coding: utf-8


import pyProBound_operator as pbo
import os
import sys

import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

from collections import Counter




print('#################################################################')
print('#### Motif generation with Probound ©Antoni Gralak 20.08.2024####')
print('#################################################################')
print('Setting env...')
this_path = os.getcwd()


# Load metadata
metadata = pd.read_csv('../exemplary_data/metadata.csv', sep=';')

data_path = '../output/00_read_in_data/'

save_path = '../output/04_pyProBound_analysis/'

try:
    os.mkdir(save_path)
except FileExistsError:
    pass



# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in sequences generated in 00_read_in and calculates a PSAM motif using ProBound.""")

parser.add_argument('-sms', '--sms_name', type=str, help='Smile-seq experiment number. E.g. exp1.', required=True)
parser.add_argument('-bs', '--binding_size', type=int, nargs='+', help="""The binding size which is to be computed by ProBound.
By default [6, 9, 12, 15, 24]. Parse multiple or a single integer.""")
parser.add_argument('-tf', '--Transcription_factor', type=str, nargs='+', help="""TFs to be included. By default all TFs
 that were approved in the experiment. Parse multiple or single TFs.""")


# Parse the command line arguments
args = parser.parse_args()
arguments = vars(args)


# Please define experiment name
experiment_name = arguments['sms_name'] 


if experiment_name in ['exp1', 'exp2', 'exp3']:
    input_id = 'input1' 
    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"
    #flanking regions of the library and corresponding parameter for ProBound model
    left = ''
    right = ''
    binding_mode_flank=5
elif experiment_name in ['exp4', 'exp5', 'exp6', 'exp7', 'exp8']:
    input_id = 'input2'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
    left = 'GGGGTACTGTGGAGATAG'
    right = 'AAACTCCCTGAGACC'
    binding_mode_flank=18
elif experiment_name in ['exp9', 'exp10', 'exp11', 'exp12', 'exp13', 'exp14']:
    input_id = 'input3'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
    left = 'GGGGTACTGTGGAGATAG'
    right = 'AAACTCCCTGAGACC'
    binding_mode_flank=18
elif experiment_name in ['exp15', 'exp16', 'exp17', 'exp18', 'exp19', 'exp20', 'exp21', 'exp22', 'exp23']:
    input_id = 'input4'
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
    left = 'GGGGTACTGTGGAGATAG'
    right = 'AAACTCCCTGAGACC'
    binding_mode_flank=18
else:
    print("-sms needs to be a experiment ID, e.g. exp1 (possible options 1 to 23). Stopping script.")
    sys.exit(1)



# For Probound, one needs to define the flanking regions. 
# I.e. left everything, the protein had contact with before the mBC, BC and lfl; and right everything after rfl. 
# Rfl is currently 5nts.

# Define what will be analyzed

if arguments['Transcription_factor']:
    to_be_analyzed = arguments['Transcription_factor']
    
else:
    to_be_analyzed = list(metadata[(metadata['experiment'] == experiment_name) & (metadata['approved'] == True)]['TF'])


if arguments['binding_size']:
    binding_mode_size = arguments['binding_size']
else:
    binding_mode_size = [6, 9, 12, 15, 24]


print('Parameters for analysis:')
print(f'Barcodes: methylated barcode {methylated_BC}; unmethylated barcode {unmethylated_BC}.')
print(f'Libraries to be analyzed: {to_be_analyzed}.')
print(f'With binding modes of the size(s): {binding_mode_size}.')


# # Loop over it, read in corresponding input and eluted samples which are stored as csvs, change all CpG of methylated sequences into mg and run ProBound


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
        if hamming_distance(input_str, s) <= 2:
            similar_strings.append(True)
        else:
            similar_strings.append(False)
    return similar_strings



print('Starting analysis!')
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
    
    
    ###################################################################################
    # change all CpGs of methylated sequences into mg
    
    # input
    mask = input_df['methl'] == methylated_BC
    
    changed_values = input_df[input_df['methl'] == methylated_BC]['random24'].str.replace('CG', 'mg', regex=False)
    
    input_df.loc[mask, 'random24'] = changed_values
    
    # eluted
    mask = eluted_df['methl'] == methylated_BC
    
    changed_values = eluted_df[eluted_df['methl'] == methylated_BC]['random24'].str.replace('CG', 'mg', regex=False)
    
    eluted_df.loc[mask, 'random24'] = changed_values
    
    ###################################################################################
    # run ProBound
    for b_size in binding_mode_size:
        
        print(f"Starting analysis for {TF} with a binding mode size of {b_size}.")
        
        # Create necessary subfolders
        # results
        try:
            os.makedirs(save_path + f'{experiment_name}/results/{TF}/')
        except FileExistsError:
            pass
        
        #try:
        #    os.mkdir(save_path + TF + f'/binding_mode_size_{b_size}/')
        #except FileExistsError:
        #    pass
        
        # misc
        try:
            os.makedirs(save_path + f'{experiment_name}/misc/{TF}/')
        except FileExistsError:
            pass

        #try:
        #    os.mkdir(save_path + TF + f'/binding_mode_size_{b_size}/')
        #except FileExistsError:
        #    pass
        
        
        results_path = save_path + f'{experiment_name}/results/{TF}/'
        misc_path = save_path + f'{experiment_name}/misc/{TF}/'

       
        barcode_flank = eluted_df[['methl', 'BC', 'lfl']].apply(lambda row: ''.join(row.values.astype(str)), axis=1)[0]
        left_flank = left + barcode_flank

        right_flank = eluted_df['rfl'][0] + right

        #creating a df for pyProBound
        df_input = pd.DataFrame(
                                {'header': np.repeat('input', len(input_df))
                                })
        df_input['sequence'] = list(input_df['random24'])

        df_eluted = pd.DataFrame(
                                {'header': np.repeat('eluted', len(eluted_df))
                                })
        df_eluted['sequence'] = list(eluted_df['random24'])


        #technically not necessary to save it with the name, this is why you have the config.alter_output func later
        outputfile = misc_path + "f_" + TF + "_output.tsv"
        count_table = pbo.build_count_table(df_input, df_eluted,
                                        output_filename=outputfile, gzip=False)

        #######################################################################
        #Config file preparation
        config = pbo.generate_meSMiLE_seq_configuration(outputfile,
                                                        24,
                                                        meCpG_encoding="mg", # this is default
                                                        left_flank=left_flank,
                                                        right_flank=right_flank,
                                                        binding_mode_flank=binding_mode_flank,
                                                        binding_modes=3,
                                                        binding_mode_size=b_size
                                                       )
        basename = TF + '_testmodel'
        config.alter_output(output_path=misc_path, 
                            base_name=basename, 
                            print_trajectory=True, 
                            verbose=False)
        config_filename = misc_path + TF + '_config.json'
        config.print_json(config_filename)

        #########################################################################
        #Run Probound and get PSAMs
        print("Running ProBound...")
        pbo.run_probound(config_filename, cleanup_verbose=False)
        
        # cleanup
        os.remove("tmp.optimization.out")

        psams = pbo.get_psam(misc_path + f"{basename}.models.json")

        #########################################################################
        #Saving and visualizing PSAMs
        for j, psam in enumerate(psams):
            psam.to_csv(results_path + TF + f"_bindingmode_{str(j + 1)}.csv")

            fig, ax = plt.subplots(1,1,figsize=[10,6])
            logo = logomaker.Logo(psams[j],
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

            fig.suptitle(f"{TF}_bindingmode_{str(j + 1)}")

            fig.savefig(results_path + TF + f"_bindingmode_{str(j + 1)}_logo.pdf", format='pdf')
            fig.savefig(results_path + TF + f"_bindingmode_{str(j + 1)}_logo.png", format='png')
            plt.close()



        print("Finished!")
print('Script finished!')


