#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyProBound_operator as pbo
import os
import sys

import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker




# In[2]:


print('#################################################################')
print('#### Motif generation with Probound ©Antoni Gralak_19.05.2025####')
print('#################################################################')
print('Setting env...')
this_path = os.getcwd()
sys.path.append(this_path)



# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in sequences generated in 00_read_in and calculates a dGG motif using ProBound.""")

parser.add_argument('-sms', '--sms_name', type=str, help='Smile-seq experiment number. E.g. SmSAG01.', required=True)
parser.add_argument('-bs', '--binding_size', type=int, nargs='+', help="""The binding size which is to be computed by ProBound.
By default [6, 9, 12, 15, 24]. Parse multiple or a single integer.""")
parser.add_argument('-bc', '--Barcode', type=str, nargs='+', help="""Barcodes to be included. By default all BC1 to 12.
Parse multiple or single BCs.""")
#parser.add_argument('-gr', '--graphs', type=str, help='Do you want summary graphs? True/False', required=True)

# Parse the command line arguments
args = parser.parse_args()
arguments = vars(args)


# In[4]:


# Please define experiment name
experiment_name = arguments['sms_name'] 


# In[5]:


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


# In[6]:


data_path = '/home/gralak/updepla/users/gralak/NAS2/SmileSeq_paper/SmileSeq_experiments/' + experiment_name + '/00_read_in_data/output/'

calc_path = f'/home/gralak/updepla/users/gralak/SmileSeq_paper/meSMiLEseq_joint_analysis/{experiment_name}/04_ProBound_analysis/calc/'
psam_path = f'/home/gralak/updepla/users/gralak/SmileSeq_paper/meSMiLEseq_joint_analysis/{experiment_name}/04_ProBound_analysis/psam/'



try:
    os.makedirs(calc_path)
except FileExistsError:
    pass

try:
    os.makedirs(psam_path)
except FileExistsError:
    pass



# In[ ]:


if experiment_name == 'SmSAG01':
    barcodes = pd.read_csv(this_path + '/' + experiment_name + '/Barcodes.csv', sep='\t', header=None)
    TFs = pd.read_csv(this_path + '/' + experiment_name + '/analysed_TFs.csv', sep=';')
    
    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"
    #flanking regions of the library and corresponding parameter for ProBound model
    left = ''
    right = ''
    binding_mode_flank=5

elif experiment_name in ['SmSAG02', 'SmSAG03']:
    barcodes = pd.read_csv(this_path + '/' + experiment_name + '/Barcodes.csv', sep=',', header=None)
    TFs = pd.read_csv(this_path + '/' + experiment_name + '/analysed_TFs.csv', sep=';')

    methylated_BC = "AGTA"
    unmethylated_BC = "GAGT"
    left = ''
    right = ''
    binding_mode_flank=5

else:
    barcodes = pd.read_csv(this_path + '/' + experiment_name + '/Barcodes.csv', header=None)
    TFs = pd.read_csv(this_path + '/' + experiment_name + '/analysed_TFs.csv')
    
    methylated_BC = "AGTA"
    unmethylated_BC = "GAAT"
    left = 'GGGGTACTGTGGAGATAG'
    right = 'AAACTCCCTGAGACC'
    binding_mode_flank=18


# ### For Probound, I need to define the flanking regions. I.e. left everything, the protein had contact with before the mBC, BC and lfl; and right everything after rfl. Rfl is currently 5nts.

# ### Barcodes and analysed TFs.

# In[8]:




TFs = TFs[barcodes.loc[:,2]]

#reindex
TFs['index'] = np.arange(0, len(TFs), step=1)
TFs = TFs.rename(index=TFs['index'])


# # Define what will be analyzed

# In[9]:


if arguments['Barcode']:
    to_be_analyzed = arguments['Barcode']
    
else:
    to_be_analyzed = ['BC1','BC2','BC3','BC4','BC5','BC6','BC7','BC8','BC9','BC10','BC11','BC12']



# In[10]:


if arguments['binding_size']:
    binding_mode_size = arguments['binding_size']
else:
    binding_mode_size = [6, 9, 12, 15, 24]


# In[ ]:


print('Parameters for analysis:')
print(f'Barcodes: methylated barcode {methylated_BC}; unmethylated barcode {unmethylated_BC}.')
print(f'Libraries to be analyzed: {to_be_analyzed}.')
print(f'With binding modes of the size(s): {binding_mode_size}.')


# # Loop over it, read in corresponding input and eluted samples which are stored as csvs, change all CpG of methylated sequences into mg and run ProBound

# In[11]:


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


# In[ ]:


print('Starting analysis!')
for BC in to_be_analyzed:
    
    print('loading necessary datasets...')
    print(f'loading {BC}...')
    input_df = pd.read_csv(input_data_path + f'{BC}_contamination_filtered.csv')
    eluted_df = pd.read_csv(data_path + f'{BC}_contamination_filtered.csv')
    
    
    
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
        name = TFs[barcodes[0] == BC]['Proteins'].values[0]
        print(f"Starting analysis for {name} with a binding mode size of {b_size}.")
        
        # Create necessary subfolders
        # results
        try:
            os.makedirs(f'{psam_path}/{name}/binding_mode_size_{b_size}/')
        except FileExistsError:
            pass
        
        #output 
        
        try:
            os.makedirs(f'{calc_path}/{name}/binding_mode_size_{b_size}/')
        except FileExistsError:
            pass

        
        psam_path_TF = f'{psam_path}/{name}/binding_mode_size_{b_size}/'
        calc_path_TF = f'{calc_path}/{name}/binding_mode_size_{b_size}/'

        # for the left flank, you of course have two different due to the mBC.. I will use mBC=AGTA for now
        # not very elegant, basically I concat the first three columns, save it as a list to be able to extract 
        # the 0th element (otherwise I need to know the index))
       
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
        outputfile = calc_path_TF + "f_" + name + "_output.tsv"  
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
        basename = name + '_testmodel'
        config.alter_output(output_path=calc_path_TF, 
                            base_name=basename, 
                            print_trajectory=True, 
                            verbose=False)
        config_filename = calc_path_TF + name + '_config.json'  
        config.print_json(config_filename)

        #########################################################################
        #Run Probound and get PSAMs
        print("Running ProBound...")
        pbo.run_probound(config_filename, cleanup_verbose=False)
        
        # ...aaaaaaand cleanup
        os.remove("tmp.optimization.out")

        psams = pbo.get_psam(calc_path_TF + f"{basename}.models.json") 

        #########################################################################
        #Saving and visualizing PSAMs
        for j, psam in enumerate(psams):
            psam.to_csv(psam_path_TF + name + f"_bindingmode_{str(j + 1)}.csv") 

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

            fig.suptitle(f"{name}_bindingmode_{str(j + 1)}")

            fig.savefig(psam_path_TF + name + f"_bindingmode_{str(j + 1)}_logo.pdf", format='pdf')
            fig.savefig(psam_path_TF + name + f"_bindingmode_{str(j + 1)}_logo.png", format='png')
            plt.close()



        print("Finished!")
print('Script finished!')






