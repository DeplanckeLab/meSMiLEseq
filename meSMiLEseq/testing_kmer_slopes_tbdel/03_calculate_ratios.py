#!/usr/bin/env python
# coding: utf-8

import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os




print('####################################################################')
print('#### Calculate ratios for scatterplots ©Antoni Gralak_20.08.2024 ###')
print('####################################################################')
print('Setting env...')
#this_path = os.getcwd()
#sys.path.append(this_path)
#sys.path.append('/data/gralak/meSMiLEseq_github/meSMiLEseq/separate_analysis/')

#import functions

# Load metadata
#metadata = pd.read_csv('../exemplary_data/metadata.csv', sep=';')

data_path = '/home/gralak/updepla/users/gralak/SmileSeq_paper/meSMiLEseq_joint_analysis/'

#save_path = '../output/03_calculate_ratios/'

#try:
#    os.mkdir(save_path)
#except FileExistsError:
#    pass


# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in sequences generated in 01_kmer_analysis and calculates kmer ratios
to generate scatterplots as used in Figure 2 in Gralak et al.""")

parser.add_argument('-sms', '--sms_name', type=str, help='Smile-seq experiment number. E.g. exp1.', required=True)
parser.add_argument('-k', '--kmer', type=int, nargs='+', help="""kmer size. By default k=[6, 7, 8, 9]. Parse multiple
 or a single integer.""")
parser.add_argument('-tf', '--Transcription_factor', type=str, nargs='+', help="""TFs to be included. By default all TFs
 that were approved in the experiment. Parse multiple or single TFs.""")

# Parse the command line arguments
args = parser.parse_args()
arguments = vars(args)

experiment_name = arguments['sms_name']

#if experiment_name in ['exp1', 'exp2', 'exp3','exp4', 'exp5', 'exp6', 'exp7', 'exp8',
#                        'exp9', 'exp10', 'exp11', 'exp12', 'exp13', 'exp14',
#                        'exp15', 'exp16', 'exp17', 'exp18', 'exp19', 'exp20', 'exp21', 'exp22', 'exp23']:
#    pass
#else:
#    print("-sms needs to be a experiment ID, e.g. exp1 (possible options 1 to 23). Stopping script.")
#    sys.exit(1)

kmer_path = os.path.join(data_path, experiment_name + '/01_kmer_analysis/')
to_be_analyzed = os.listdir(kmer_path)

save_path = os.path.join(data_path, experiment_name + '/03_kmer_ratios/')

try:
    os.makedirs(save_path)
except FileExistsError:
    pass

# Define what will be analyzed

#if arguments['Transcription_factor']:
#    to_be_analyzed = arguments['Transcription_factor']
    
#else:
#    to_be_analyzed = list(metadata[(metadata['experiment'] == experiment_name) & (metadata['approved'] == True)]['TF'])


if arguments['kmer']:
    kmers = arguments['kmer']
else:
    kmers = [6, 7, 8, 9]


for file in to_be_analyzed:
    TF = file.split('_')[0] + ' ' + file.split('_')[1]
    for kmer in kmers:
        print(f'Calculating ratios for {TF} ({experiment_name}, {kmer}mer) with scatterplot..')
        #ratio_df = functions.calculate_ratios(experiment_name, TF, kmer, data_path)

        kmer_df = pd.read_csv(os.path.join(kmer_path, file))
        kmer_e = kmer_df[kmer_df['status'] == 'eluted']
        pattern = '|'.join(kmer_e['kmer'])
        filtered_kmer = kmer_df[kmer_df['kmer'].str.contains(pattern, regex=True)]

        m_dict = {}
        nm_dict = {}
        for k, df in filtered_kmer.groupby('kmer'):
            # extract count of methylated kmer in input
            im_count = df[(df['mod'] == 'methl') & (df['status'] == 'input')]['count'].values
            # generate ratio and save as dictionary
            m_dict[k] = df[df['mod'] == 'methl']['count'].values/im_count
            # extract count of nonmethylated kmer in input
            inm_count = df[(df['mod'] == 'nonmethl') & (df['status'] == 'input')]['count'].values
            # generate ratio
            nm_dict[k] = df[df['mod'] == 'nonmethl']['count'].values/inm_count

        
    
        new_df = pd.DataFrame.from_dict(m_dict, orient='index', columns=['input_methl', 'eluted_methl']).T.replace([], np.nan)
        new_df = new_df.T.reset_index()

        new_df2 = pd.DataFrame.from_dict(nm_dict, orient='index', columns=['input_nonmethl', 'eluted_nonmethl']).T.replace([], np.nan)
        new_df2 = new_df2.T.reset_index()
            
            #add CG information
        final_df = pd.merge(new_df, new_df2, on='index')
        final_df['CpG'] = ['CG' in final_df['index'][i] for i in range(len(final_df))]

        name_file = file.replace('enrichment', 'ratio')
        final_df.to_csv(os.path.join(save_path, name_file), index=False)

        
        # Plot scatterplot using ratios
        
        # Some colors to choose from
        black = '#000000'
        lightblack = '#333333'
        darkgray = '#666666'
        mediumgray = '#999999'
        lightgray = '#CCCCCC'

        darkred = '#FF0000'
        red = '#FF3333'
        lightred = '#FF6666'
        pink ='#FF9999'
        salmon = '#FFCCCC'



        x = final_df['eluted_methl']
        y = final_df['eluted_nonmethl']
        
        
        fig, ax = plt.subplots(1, 1, figsize=(5.5, 5.5))

        #colors = {'present' : darkred, 'not present' : lightblack}
        colors = {True : darkred, False : lightblack}

        ax.scatter(x=x, y=y, alpha=0.6, facecolors='none', edgecolors=final_df['CpG'].map(colors), rasterized = True)
        ax.grid(visible=False)


        # Remove the top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        lower_limit = None
        # Extent the axis by 5 % of max value
        
        if x.nlargest(1).values[0] > y.nlargest(1).values[0]:
            ax.set_xlim(lower_limit,x.nlargest(1).values[0] + (0.05*x.nlargest(1).values[0]))
            ax.set_ylim(lower_limit,x.nlargest(1).values[0] + (0.05*x.nlargest(1).values[0]))

        else:
            ax.set_xlim(lower_limit,y.nlargest(1).values[0]+(0.05*y.nlargest(1).values[0]))
            ax.set_ylim(lower_limit,y.nlargest(1).values[0]+(0.05*y.nlargest(1).values[0]))

        ax.set_aspect('equal', adjustable='box')
        
        
        ax.set_xlabel('methylated DNA', fontfamily='sans-serif', fontsize=10, fontstyle='italic')
        ax.set_ylabel('unmethylated DNA', fontfamily='sans-serif', fontsize=10, fontstyle='italic')
        
        ax.set_title(f"{experiment_name}, {TF}, {kmer} enrichment, normalized by input")

        #ax.ticklabel_format(style='sci', scilimits=(0,0))

        plt.savefig(save_path + f'/{TF}_{kmer}mer_scatterplot.pdf')

        print('Done!')






