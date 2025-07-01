#!/usr/bin/env python
# coding: utf-8

# In[170]:


import argparse
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.lines
import numpy as np
import csv


# In[ ]:


print('#############################################################')
print('#### Creating kmer scatterplots ©Antoni Gralak_30.05.2023####')
print('#############################################################')
print('Setting env...')
this_path = os.getcwd()
sys.path.append(this_path)


# In[ ]:


# Creating a parser argument
parser = argparse.ArgumentParser("""This script reads in data generated in 01_kmer_analysis and creates scatterplots.""")

parser.add_argument('-sms', '--sms_name', type=str, help='Smile-seq experiment number. E.g. SmSAG01.', required=True)
parser.add_argument('-k', '--kmer', type=int, help="""kmer size. By default k=8. Parse a single integer.""")
parser.add_argument('-tf', '--TF', type=str, nargs='+', help="""TFs to be included. By default all TFs from the experiment.
#Parse multiple or single BCs.""")
parser.add_argument('-p', '--p_value', type=float, help='Set a p-value as threshold. By default 0.01')
parser.add_argument('-f', '--filter', type=int, help='Do you want to filter by p value. By default 1 (= True).')
parser.add_argument('-s', '--slope', type=int, help='Plot slope. By default 1 (= True).')
parser.add_argument('-log', '--log', type=int, help='If log then log2. By default 0 (= False).')

# Parse the command line arguments
args = parser.parse_args()
arguments = vars(args)


# In[ ]:


experiment_name = arguments['sms_name'] 


# In[ ]:


#setting input to correct for library bias.
if experiment_name in ['SmSAG01', 'SmSAG02', 'SmSAG03']:
    input_data_path = '/home/gralak/NAS2/gralak/SmileSeq_paper/SmileSeq_experiments/inputs/20220315_input/00_read_in_data/output/'
elif experiment_name in ['SmSAG04', 'SmSAG05', 'SmSAG06', 'SmSAG07', 'SmSAG08']:
    input_data_path = '/home/gralak/NAS2/gralak/SmileSeq_paper/SmileSeq_experiments/inputs/20220915_input/00_read_in_data/output/'
elif experiment_name in ['SmSAG09', 'SmSAG10', 'SmSAG11', 'SmSAG12', 'SmSAG13', 'SmSAG14']:
    input_data_path = '/home/gralak/NAS2/gralak/SmileSeq_paper/SmileSeq_experiments/inputs/20230228_input/00_read_in_data/output/'
elif experiment_name in ['SmSAG15', 'SmSAG16', 'SmSAG17', 'SmSAG18', 'SmSAG19', 'SmSAG20', 'SmSAG21', 'SmSAG22', 'SmSAG23']:
    input_data_path = '/home/gralak/NAS2/gralak/SmileSeq_paper/SmileSeq_experiments/inputs/20230503_input/00_read_in_data/output/'
else:
    print("-sms needs to be a SmileSeq name such as SmSAG01 (possible options 01 to 23). Stopping script.")
    sys.exit(1)


# In[ ]:


# setting paths to obtain and save data (data_ save_path, respecitvely)
data_path = '/home/gralak/NAS2/gralak/SmileSeq_paper/SmileSeq_experiments/' + experiment_name + '/01_kmer_analysis/output/'

save_path = '/home/gralak/NAS2/gralak/SmileSeq_paper/plotting/kmer_scatterplots/'


# In[ ]:


#getting info about abnalysed proteins and setting barcodes appropriatly
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




# In[ ]:


#setting parameters, i.e. Which TFs to analyse and which kmers to use and p-values
if arguments['TF']:
    to_be_analyzed = arguments['TF']
    
else:
    to_be_analyzed = list(TFs['Proteins'])


if arguments['kmer']:
    kmer = arguments['kmer']
else:
    kmer = 8

if arguments['p_value']:
    p_val = arguments['p_value']
else:
    p_val = 0.01

if arguments['filter'] == 0:
    fltr = arguments['filter']
else:
    fltr = 1

if arguments['slope'] == 0:
    slp = arguments['slope']
else:
    slp = 1

if arguments['log'] == 1:
    lg = 1
else:
    lg = 0


# In[ ]:


# Used colorpalette
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


# In[270]:


for iterator, TF in enumerate(to_be_analyzed):
    #Create folders to save everything
    try:
        os.mkdir(save_path + f'{TF}/')
    except FileExistsError:
        pass

    output_path = save_path + f'{TF}/'

    print('Loading datasets')
    print(TF)
    
    BC = TFs[TFs['Proteins'] == TF]
    BC = BC['Mitomi Buttons'].values[0]

    kmer_df = pd.read_csv(data_path + f'{TF}_{kmer}mer_enrichment.csv')
    pval_df = pd.read_csv(data_path + f'{TF}_{kmer}_pvalues.csv')
    
    # loading corresponding input and calculating bias within the library
    i_df = pd.read_csv(input_data_path + f'BC{BC+1}_contamination_filtered.csv')
    total = len(i_df)
    m = len(i_df[i_df['methl'] == methylated_BC])
    nm = len(i_df[i_df['methl'] == unmethylated_BC])

    m_bias = m/total
    nm_bias = nm/total
    
    # filter p_val for significant values
    pval_df = pval_df[pval_df['p_adjust'] < p_val]
    
    # Divide and conquer, split by methylation and filter the kmer_df
    p_val_m = pval_df[pval_df['mod'] == 'methl']
    p_val_nm = pval_df[pval_df['mod'] == 'nonmethl']
    
    #If filter is true then filtering by pvalue, if not, not.
    if fltr == 1:
        methl = kmer_df[(kmer_df['mod'] == 'methl') & (kmer_df['status'] == 'eluted')]
        methl = methl[methl['kmer'].isin(p_val_m['kmer'])]
        
        nonmethl = kmer_df[(kmer_df['mod'] == 'nonmethl') & (kmer_df['status'] == 'eluted')]
        nonmethl = nonmethl[nonmethl['kmer'].isin(p_val_nm['kmer'])]

    elif fltr == 0:
        methl = kmer_df[(kmer_df['mod'] == 'methl') & (kmer_df['status'] == 'eluted')]
        nonmethl = kmer_df[(kmer_df['mod'] == 'nonmethl') & (kmer_df['status'] == 'eluted')]
    
    # sanity check
    if (len(methl) < 5) | (len(nonmethl) < 5):
        print('Not enough reads to plot scatterplot.')
        print(f'Methylated_df: {len(methl)}')
        print(f'Unmethylated_df: {len(nonmethl)}')
        continue
    
    
    # adjust counts for input_bias
    methl['count'] = methl['count']/m_bias
    nonmethl['count'] = nonmethl['count']/nm_bias
    
    
    
    print('Plotting!')
    
    # Plot everything
    CpG_m = methl[['kmer', 'CpG']]
    CpG_nm = nonmethl[['kmer', 'CpG']]
    CpG = pd.concat([CpG_m, CpG_nm])

    methl = methl[['kmer', 'count']]
    nonmethl = nonmethl[['kmer', 'count']]

    combine = pd.merge(methl, nonmethl, how='outer', on='kmer')
    combine = combine.fillna(0)


    combine = pd.merge(combine, CpG, how='outer', on='kmer')

    if lg == 1:
        combine['count_x'] = combine['count_x'].apply(lambda x: np.log2(x))
        combine['count_y'] = combine['count_y'].apply(lambda x: np.log2(x))
        
        combine = combine[(combine['count_x'] > 5 ) & (combine['count_y'] > 5)]

    

    fit_CpG = combine[combine['CpG'] == 'present']
    fit_no_CpG = combine[combine['CpG'] == 'not present']
    
    #for linear regression calculation exclude all kmers which are zero in either x or y value
    fit_CpG = fit_CpG[(fit_CpG['count_x'] != 0) & (fit_CpG['count_y'] != 0)]
    fit_no_CpG = fit_no_CpG[(fit_no_CpG['count_x'] != 0) & (fit_no_CpG['count_y'] != 0)]
  

    x = combine['count_x']
    y = combine['count_y']
    
   

    # create scatterplot

    fig, ax = plt.subplots(1, 1, figsize=(5.5, 5.5))


    colors = {'present' : darkred, 'not present' : lightblack}

    ax.scatter(x=x, y=y, alpha=0.6, facecolors='none', edgecolors=combine['CpG'].map(colors))
    ax.grid(visible=False)




    # Remove the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    if slp == 1:

        # Linear regression line for methylated reads
        if len(fit_CpG) > 0:
            slope, intercept = np.polyfit(fit_CpG['count_x'], fit_CpG['count_y'], 1)
            x_values = np.linspace(min(fit_CpG['count_x']), max(fit_CpG['count_x']))
            y_values = slope * x_values + intercept
            ax.plot(x_values, y_values, color=darkred, linewidth=1)
            # Create the equation string
            equation = f"f(x) = {slope:.2f}x + {intercept:.2f}"
        else:
            equation = slope = 'no kmers with CpG'

        
        
        # Linear regression line for unmethylated reads
        if len(fit_no_CpG) > 0:
            slope2, intercept2 = np.polyfit(fit_no_CpG['count_x'], fit_no_CpG['count_y'], 1)
            x_values2 = np.linspace(min(fit_no_CpG['count_x']), max(fit_no_CpG['count_x']))
            y_values2 = slope2 * x_values2 + intercept2
            ax.plot(x_values2, y_values2, color=black, linewidth=1)
            # Create the equation string
            equation2 = f"f(x) = {slope2:.2f}x + {intercept2:.2f}"
        else:
            equation2 = slope2 = 'no kmers without CpG'
    else:
        slope = slope2 ='not available'

    

    if lg == 0:
        lower_limit = None
    else:
        if x.nsmallest(1).values[0] < y.nsmallest(1).values[0]:
            lower_limit = x.nsmallest(1).values[0] - 0.5
        else:
            lower_limit = y.nsmallest(1).values[0] - 0.5

    # Extent the axis by 5 % of max value
    if x.nlargest(1).values[0] > y.nlargest(1).values[0]:
        ax.set_xlim(lower_limit,x.nlargest(1).values[0] + (0.05*x.nlargest(1).values[0]))
        ax.set_ylim(lower_limit,x.nlargest(1).values[0] + (0.05*x.nlargest(1).values[0]))

    else:
        ax.set_xlim(lower_limit,y.nlargest(1).values[0]+(0.05*y.nlargest(1).values[0]))
        ax.set_ylim(lower_limit,y.nlargest(1).values[0]+(0.05*y.nlargest(1).values[0]))

    ax.set_aspect('equal', adjustable='box')    

    if slp == 1:
        ax.set_xlabel(f'methylated DNA\n red: {equation}, black: {equation2}', 
                    fontfamily='sans-serif', fontsize=10, fontstyle='italic')
    else:
        ax.set_xlabel(f'methylated DNA', 
                    fontfamily='sans-serif', fontsize=10, fontstyle='italic')
    ax.set_ylabel('unmethylated DNA', fontfamily='sans-serif', fontsize=10, fontstyle='italic')

    if lg == 0:
        ax.ticklabel_format(style='sci', scilimits=(0,0))
    
    #Setting title 
    plt.title(f'{TF}: {experiment_name}')

    # Set rasterization to True
    ax.set_rasterized(True)

    plt.savefig(output_path + f'{experiment_name}_{kmer}mer.png', transparent=True, dpi=400)
    plt.savefig(output_path + f'{experiment_name}_{kmer}mer.pdf', transparent=True, dpi=400)
    plt.savefig(output_path + f'{experiment_name}_{kmer}mer.svg', transparent=True, dpi=500)

    
    #saving all parameters in csv
    print('Saving slopes, most significant kmers etc..')
    
    
    file_path = os.path.join(output_path, f'{kmer}_specifications.csv')
    
    try:
        most_significant = pval_df.nsmallest(columns='p_adjust', n=1)['kmer'].values[0]
        pvalue = pval_df.nsmallest(columns='p_adjust', n=1)['p_adjust'].values[0]
        mod = pval_df.nsmallest(columns='p_adjust', n=1)['mod'].values[0]
    except IndexError:
        most_significant = pvalue = mod = 'not available'
  

    
    
    
    specs = [experiment_name, most_significant, pvalue, mod, slope, slope2]


    # Check if the file exists
    if os.path.exists(file_path):
        # Append data to the existing file
        with open(file_path, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(specs)
    else:
        # Create a new file and write data
        with open(file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['SmS', 'most_significant', 'p_val', 'mod', 'CpG_slope', 'no_CpG_slope'])
            writer.writerow(specs)
    
    print('Finished!')







