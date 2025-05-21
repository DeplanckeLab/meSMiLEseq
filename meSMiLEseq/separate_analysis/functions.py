#############################################################################
#### Functions for Data Wrangling with SmileSeq sequences ©Antoni Gralak ####
#############################################################################

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.lines
import collections
import os

import itertools

######################################################################################
# Most important function is kmer_counting. All other functions                      #
# are somewhat related to this one. kmer_counting takes sequencing data from one or  #
# multiple DFs with sequences (sequences are rows, columns must include one column   # 
# called 'methl') and divides the sequence into predefined kmers. It generates a tidy#
# dataframe with columns: 'kmer', 'count' (being trivial), 'status' -> 'input', 'FT' #
# 'eluted' (or something completely different like BC-specific, green, cows), 'mod'  #
# -> 'mod' is meant to represent the modification of the DNA. Gralak et al worked    #
# with methylated and nonmethylated DNA, so for now mod is 'methl' or 'nonmethl'.    #
#                                                                                    #
#                                                                                    #
######################################################################################



def filter_my_files_for_fastq(all_files, pattern):
    """This function saves all filenames in a list which have the defined pattern in their name."""
    f1 = str(pattern)
    f2 = re.compile(f1+'(.*?)\.fastq')
    ll = []
    for x in all_files:
        found = re.search(f2, x)
        if found:
            ll.append(x)
    return(ll)




def readin_fastq(core_path,filename):
    """This function reads in fastq files and saves them as pd.df."""
    seqsraw = pd.read_csv(os.path.join(core_path,filename),sep="\t",header=None)  # txt files are tab-separated that is why we say that the seperator is "tab" represented by "\t" (sep = "\t")
    seqs = seqsraw.iloc[range(1,seqsraw.shape[0],4),:]
    return(seqs)




def rev_complement(seq):
    """Self-explanatory, creates reverse complement."""
    thisdict = {
    "G": "C",
    "C": "G",
    "A": "T",
    "T": "A",
    "N": "N",
    "M": "G"
    }
    translated = ''
    for k in seq:
        translated = translated + thisdict[k]
    reverse_translated = translated[::-1]
    return reverse_translated   


# kmer counting
#This function divides the 24 random region into 8mers (per definit) and counts (uses the function Counter) how often each kmer appears -> helps to identify bias in the library itself

def kmer_counting(df, kmer = 8, status = None, one_df = True, mbc = ['AGTA','GAAT']):
    """This function divides 24 random nucleotides from methlSmileSeq libraries into kmers of a specified, 
    changable length. 
    Return a tidy df showing status = [elute, input]; mod = [methl, nonmethl]; kmer = extracted kmer; 
    count = count of kmer.
    
    You can either pass this function one single df with only one status, then set one_df = True and define status as
    character string. 
    Or you can pass it a list of dfs with different conditions (e.g. input, elution, flow through). 
    Then set one_df = False and pass status a list of conditions."""
    
    if one_df == True:
        methl = df[df['methl'] == mbc[0]] #split for methylated BC
        methl.index = range(0,len(methl))
        nonmethl = df[df['methl'] == mbc[1]]
        nonmethl.index = range(0,len(nonmethl))
    
        methl_kmers = [methl['random24'][i][x:(x+kmer)] for i in range(len(methl)) for x in range(24-kmer+1)] #generates all kmers
        nonmethl_kmers = [nonmethl['random24'][i][x:(x+kmer)] for i in range(len(nonmethl)) for x in range(24-kmer+1)]
    
        methl_kmers = pd.DataFrame.from_dict(dict(collections.Counter(methl_kmers)),orient='Index') #counts kmers
        nonmethl_kmers = pd.DataFrame.from_dict(dict(collections.Counter(nonmethl_kmers)),orient='Index')
        
        methl_kmers = methl_kmers.reset_index()
        nonmethl_kmers = nonmethl_kmers.reset_index()

        methl_kmers = methl_kmers.rename(columns={0:'count', 'index':'kmer'})
        nonmethl_kmers = nonmethl_kmers.rename(columns={0:'count', 'index':'kmer'})

        methl_kmers['mod'] = 'methl'
        nonmethl_kmers['mod'] = 'nonmethl'

        methl_kmers['status'] = status
        nonmethl_kmers['status'] = status

        kmer_df = pd.concat([methl_kmers, nonmethl_kmers])
        kmer_df = kmer_df.reset_index()
        kmer_df = kmer_df.drop(columns=['index'])
        return kmer_df
    
    elif one_df == False:
        kmer_df = pd.DataFrame()
        for i,my_df in enumerate(df):
            methl = my_df[my_df['methl'] == mbc[0]] #split for methylated BC
            methl.index = range(0,len(methl))
            nonmethl = my_df[my_df['methl'] == mbc[1]]
            nonmethl.index = range(0,len(nonmethl))

            methl_kmers = [methl['random24'][i][x:(x+kmer)] for i in range(len(methl)) for x in range(24-kmer+1)] #generates all kmers
            nonmethl_kmers = [nonmethl['random24'][i][x:(x+kmer)] for i in range(len(nonmethl)) for x in range(24-kmer+1)]

            methl_kmers = pd.DataFrame.from_dict(dict(collections.Counter(methl_kmers)),orient='Index') #counts kmers
            nonmethl_kmers = pd.DataFrame.from_dict(dict(collections.Counter(nonmethl_kmers)),orient='Index')

            methl_kmers = methl_kmers.reset_index()
            nonmethl_kmers = nonmethl_kmers.reset_index()

            methl_kmers = methl_kmers.rename(columns={0:'count', 'index':'kmer'})
            nonmethl_kmers = nonmethl_kmers.rename(columns={0:'count', 'index':'kmer'})

            methl_kmers['mod'] = 'methl'
            nonmethl_kmers['mod'] = 'nonmethl'

            methl_kmers['status'] = status[i]
            nonmethl_kmers['status'] = status[i]

            sub_kmer_df = pd.concat([methl_kmers, nonmethl_kmers])
            sub_kmer_df = sub_kmer_df.reset_index()
            sub_kmer_df = sub_kmer_df.drop(columns=['index'])
            
            kmer_df = pd.concat([kmer_df, sub_kmer_df])
            kmer_df = kmer_df.reset_index()
            kmer_df = kmer_df.drop(columns=['index'])
        return kmer_df

# function, calculate ratios

def calculate_ratios(experiment_name, TF, kmer, data_path):
    """This function takes kmer csvs and 'normalizes' (i.e. creates ratios eluted/input) 
    the eluted kmers with the input. Required: experiment_name, TF, kmer, data_path."""
    experiment_name = experiment_name 
    TF = TF 
    kmer = kmer
    data_path = data_path
    log = 'normalised by input'

    
    kmer_df = pd.read_csv(data_path + f'{experiment_name}_{TF}_{kmer}mer_enrichment.csv')

    
    
    kmer_e = kmer_df[kmer_df['status'] == 'eluted']
    kmer_i = kmer_df[kmer_df['status'] == 'input']

    
        
    # write all kmers as a regex pattern
    pattern = '|'.join(kmer_e['kmer'])

    # Use .str.contains with regex=True to filter kmer_df for all eluted kmers
    filtered_kmer = kmer_df[kmer_df['kmer'].str.contains(pattern, regex=True)]

    
    
    # iterate over the filtered kmer dataframe and divide the eluted count by input 
    # counts (methylated and unmethylated separately)
    m_dict = {}
    nm_dict = {}
    for df in filtered_kmer.groupby('kmer'):
        # extract count of methylated kmer in input
        im_count = df[1][(df[1]['mod'] == 'methl') & (df[1]['status'] == 'input')]['count'].values
        # generate ratio and save as dictionary
        m_dict[df[0]] = df[1][df[1]['mod'] == 'methl']['count'].values/im_count
        # extract count of nonmethylated kmer in input
        inm_count = df[1][(df[1]['mod'] == 'nonmethl') & (df[1]['status'] == 'input')]['count'].values
        # generate ratio
        nm_dict[df[0]] = df[1][df[1]['mod'] == 'nonmethl']['count'].values/inm_count

    # store everything in a new df
    new_df = pd.DataFrame.from_dict(m_dict, orient='index', columns=['input_methl', 'eluted_methl']).T.replace([], np.nan)
    new_df = new_df.T.reset_index()

    new_df2 = pd.DataFrame.from_dict(nm_dict, orient='index', columns=['input_nonmethl', 'eluted_nonmethl']).T.replace([], np.nan)
    new_df2 = new_df2.T.reset_index()
    
    #add CG information
    final_df = pd.merge(new_df, new_df2, on='index')
    final_df['CpG'] = ['CG' in final_df['index'][i] for i in range(len(final_df))]
    

    return final_df


# function, create ecdf
def kmer_ecdf(df, scale='', path='', save=True, status=None, palette=['#8da0cb','#fc8d62','#66c2a5'], TF = ''):
    """This function is an immediate follow-up of the function kmer_counting. The output, a tidy dataframe containing at least following columns: 
    'kmer', 'count', 'mod', 'status'; more won't hurt. Please define status as string or as a list of strings and it will iterate over it. 

    Please define scale. Possible options are ['raw', 'log2', 'log10'].
    Please name the analysed TF.
    Please define path where figures will be saved.

    Feel free to define your colorpalette, the used standard has three neat colors.
    """
    if status == None:
        print('Please define status!')
        return

    kmer = len(df['kmer'][0])
    #create plt using matplotlib and save in results

    #iteartion over list with status, e.g. input, eluted, FT
    if isinstance(status, list) == True:

        fig, ax = plt.subplots()

        for position in range(len(status)):

            x = pd.Series(df[df['status'] == status[position]]['count'])

            if scale == 'raw':
                x = x
            elif scale == 'log2':
                x = np.log2(x)
            elif scale == 'log10':
                x = np.log10(x)

            y = x.rank(method='first') / len(x)
            ax.scatter(x=x, y=y, c=palette[position], label=status[position])
    
        ax.set_ylabel('ECDF(kmer_count)')
        ax.set_xlabel(f'{scale}_kmer_count')
        ax.legend()
        ax.grid(visible=True, alpha=0.2)

        plt.yticks(np.arange(0, 1.1, step=0.2))
        plt.title(TF)

        if save:
            fig.savefig(path + '_' + str(kmer) + 'mer_ecdf.pdf', format='pdf')
            fig.savefig(path + '_' + str(kmer) + 'mer_ecdf.tif', format='tif')
            plt.close()
            return
        else:
            plt.show()
            return


    else:
        fig, ax = plt.subplots()

        x = pd.Series(df[df['status'] == status]['count'])

        if scale == 'raw':
            x = x
        elif scale == 'log2':
            x = np.log2(x)
        elif scale == 'log10':
            x = np.log10(x)

        y = x.rank(method='first') / len(x)
        ax.scatter(x=x, y=y, c=palette[0], label=status)

        ax.set_ylabel('ECDF(kmer_count)')
        ax.set_xlabel(f'{scale}_kmer_count')
        ax.legend()
        ax.grid(visible=True, alpha=0.2)

        plt.yticks(np.arange(0, 1.1, step=0.2))
        plt.title(TF)


        if save:
            fig.savefig(path + '_' + str(kmer) + 'mer_ecdf.pdf', format='pdf')
            fig.savefig(path + '_' + str(kmer) + 'mer_ecdf.tif', format='tif')
            plt.close()
            return
        else:
            plt.show()
            return

   



# normalise by input and add CpG info


# If extracting something from a df, it is a pd.Series object -> has an index and a value. Therefore, to access the value
# type xxx.values[0] <-- returns the zeroth value of the object. Indices are useful for merging with largelist


def normalise_by_status(df, status, CpG=True):
    """This function assumes a tidy dataframe as input with following columns: 'kmer', 'count', 'mod', 'status'. Important to notice,
    that 'status' must be defined as a string in the tidy dataframe, e.g. 'input', 'flow-through', since these values will be taken to normalise
    all entries. Equally, mod must be at least ['methl' and 'nonmethl']. Also it appends information whether CG is in the kmer. To disable that function, set CpG=False."""

    DF = pd.DataFrame(columns=['normalised_counts'])

    for kmer_seq, kmer_df in df.groupby('kmer'):

        norm_count_m = kmer_df[(kmer_df['status'] == status) & (kmer_df['mod'] == 'methl')]['count']
        norm_count_nm = kmer_df[(kmer_df['status'] == status) & (kmer_df['mod'] == 'nonmethl')]['count']

        if norm_count_m.empty:
            norm_count_m = pd.Series([0])
        if norm_count_nm.empty:
            norm_count_nm = pd.Series([0])

        df1 = kmer_df[kmer_df['mod'] == 'methl']['count'].div(norm_count_m.values[0]).to_frame(name='normalised_counts')
        df2 = kmer_df[kmer_df['mod'] == 'nonmethl']['count'].div(norm_count_nm.values[0]).to_frame(name='normalised_counts')

        if CpG:
            if 'CG' in kmer_seq:
                df1['CpG'] = 'present'
                df2['CpG'] = 'present'
            else:
                df1['CpG'] = 'not present'
                df2['CpG'] = 'not present'

        DF = pd.concat([DF,df1])
        DF = pd.concat([DF,df2])
    
    return DF





# create scatterplot for SmileSeq data
def kmer_scatter(df, scale, path='', column='count', save=True, TF = '', normalised=True):
    """Self explanatory, function creates scatterplots with x axis being methylated and y axis being nonmethylated. Colorcoding via CpG information.
    Column can be either 'count' or 'normalised_counts'.
    
    Please define path to save figure.
    Please define scale. Possible options are 'raw', 'log2' or 'log10'."""

    kmer = len(df['kmer'][0])

    if normalised:


        if column == 'normalised_counts':   
            #filter out all infs
            df = df[df['normalised_counts'] != np.inf]


        methl = df[(df['mod'] == 'methl') & (df['status'] == 'eluted')]
        nonmethl = df[(df['mod'] == 'nonmethl') & (df['status'] == 'eluted')]

        # Here I extract the CpG information. I have to do it in this seperate step since otherwise 
        # if done later with 'combined', I'd have three stats: present, not present and 0. 
        
        CpG_m = methl[['kmer', 'CpG']]
        CpG_nm = nonmethl[['kmer', 'CpG']]
        CpG = pd.concat([CpG_m, CpG_nm])
        
        
        methl = methl[['kmer', 'count', 'normalised_counts']]
        nonmethl = nonmethl[['kmer', 'count', 'normalised_counts']]
        
        # combining methl and nonmethl in this order will result in methl = x axis with _x in the colname and nonmethl
        # on the y axis with _y in colnames. Also remove kmers which only appear in either methylated or nonmethylated
        # DNA (thus how='inner').
        #EDIT: I understand the thought behind my own code, that I filter out all kmers which appear only in one of the two 
        #libraries. But what if only methylated libraries are bound? Obviously there will be kmers which don't have a CpG
        # but still were part of the bound library. Thus I should allow a unique status and fill NaN with 0.


        #combine = pd.merge(methl, nonmethl, how='inner', on='kmer')
        combine = pd.merge(methl, nonmethl, how='outer', on='kmer')
        combine = combine.fillna(0)


        combine = pd.merge(combine, CpG, how='outer', on='kmer')

        combine = combine.loc[~((combine[column + '_x'] < 0.999) & (combine[column + '_y'] < 0.999)),:]

        if scale == 'raw':

            x = combine[column + '_x']
            y = combine[column + '_y']
            lower_limit = 0.1
    
        elif scale == 'log2':

            #add pseudocounts and set lower limit
            combine[column + '_x'] += 0.1
            combine[column + '_y'] += 0.1             
            x = np.log2(combine[column + '_x'])
            y = np.log2(combine[column + '_y'])
            lower_limit = -2

        elif scale == 'log10':

            #add pseudocounts and set lower limit
            combine[column + '_x'] += 0.1
            combine[column + '_y'] += 0.1       
            x = np.log10(combine[column + '_x'])
            y = np.log10(combine[column + '_y'])
            lower_limit = -1
        
        
        # create scatterplot

        fig, ax = plt.subplots()
        

        colors = {'present' : 'C1', 'not present' : 'C0'}

        ax.scatter(x=x, y=y, alpha=0.6, facecolors='none', edgecolors=combine['CpG'].map(colors))
        ax.grid(visible=True, alpha=0.2)
        
        ax.set_xlabel(f'{scale}_methylated_DNA')
        ax.set_ylabel(f'{scale}_not_methylated_DNA')
        plt.title(TF + ' ' + column)

        if (x.nlargest(1).empty) | (y.nlargest(1).empty):
            pass
        else:
            ax.set_xlim(lower_limit - 1,x.nlargest(1).values[0] + 1)
            ax.set_ylim(lower_limit - 1,y.nlargest(1).values[0] + 1)
        

        
        # create a custom legend
        custom_legend = [matplotlib.lines.Line2D([], [], marker='o', markersize=8, color='C1', markerfacecolor='None', linestyle='None'),
                        matplotlib.lines.Line2D([], [], marker='o', markersize=8, color='C0', markerfacecolor='None', linestyle='None')]
        
        ax.legend(custom_legend, ['present', 'not present'], loc='best', fontsize=12)

        if save:
            fig.savefig(path + column + '_' + str(kmer) + 'mer_scatterplot.pdf', format='pdf')
            fig.savefig(path + column + '_' + str(kmer) + 'mer_scatterplot.tif', format='tif')
            plt.close()
            return
        else:
            plt.show()    
            return
    
    else:

        methl = df[(df['mod'] == 'methl') & (df['status'] == 'eluted')]
        nonmethl = df[(df['mod'] == 'nonmethl') & (df['status'] == 'eluted')]

        # Here I extract the CpG information. I have to do it in this seperate step since otherwise 
        # if done later with 'combined', I'd have three stats: present, not present and 0. 
        
        CpG_m = methl[['kmer', 'CpG']]
        CpG_nm = nonmethl[['kmer', 'CpG']]
        CpG = pd.concat([CpG_m, CpG_nm])
        
        
        methl = methl[['kmer', 'count']]
        nonmethl = nonmethl[['kmer', 'count']]
        
        # combining methl and nonmethl in this order will result in methl = x axis with _x in the colname and nonmethl
        # on the y axis with _y in colnames. Also remove kmers which only appear in either methylated or nonmethylated
        # DNA (thus how='inner').
        #EDIT: I understand the thought behind my own code, that I filter out all kmers which appear only in one of the two 
        #libraries. But what if only methylated libraries are bound? Obviously there will be kmers which don't have a CpG
        # but still were part of the bound library. Thus I should allow a unique status and fill NaN with 0.


        #combine = pd.merge(methl, nonmethl, how='inner', on='kmer')
        combine = pd.merge(methl, nonmethl, how='outer', on='kmer')
        combine = combine.fillna(0)
    
        
        combine = pd.merge(combine, CpG, how='outer', on='kmer')


        if scale == 'raw':

            x = combine[column + '_x']
            y = combine[column + '_y']
            lower_limit = 0
    
        elif scale == 'log2':

            #add pseudocount
            combine[column + '_x'] += 0.1
            combine[column + '_y'] += 0.1        
            x = np.log2(combine[column + '_x'])
            y = np.log2(combine[column + '_y'])
            lower_limit = -2

        elif scale == 'log10':

            #add pseudocount
            combine[column + '_x'] += 0.1
            combine[column + '_y'] += 0.1        
            x = np.log10(combine[column + '_x'])
            y = np.log10(combine[column + '_y'])
            lower_limit = -1
        
        
        # create scatterplot

        fig, ax = plt.subplots()
        

        colors = {'present' : 'C1', 'not present' : 'C0'}

        ax.scatter(x=x, y=y, alpha=0.6, facecolors='none', edgecolors=combine['CpG'].map(colors))
        ax.grid(visible=True, alpha=0.2)
        
        ax.set_xlabel(f'{scale}_methylated_DNA')
        ax.set_ylabel(f'{scale}_not_methylated_DNA')
        plt.title(TF + ' ' + column)

        if (x.nlargest(1).empty) | (y.nlargest(1).empty):
            pass
        else:
            ax.set_xlim(lower_limit - 1,x.nlargest(1).values[0] + 1)
            ax.set_ylim(lower_limit - 1,y.nlargest(1).values[0] + 1)
        
        
        # create a custom legend
        custom_legend = [matplotlib.lines.Line2D([], [], marker='o', markersize=8, color='C1', markerfacecolor='None', linestyle='None'),
                        matplotlib.lines.Line2D([], [], marker='o', markersize=8, color='C0', markerfacecolor='None', linestyle='None')]
        
        ax.legend(custom_legend, ['present', 'not present'], loc='best', fontsize=12)

        if save:
            fig.savefig(path + column + '_' + str(kmer) + 'mer_scatterplot.pdf', format='pdf')
            fig.savefig(path + column + '_' + str(kmer) + 'mer_scatterplot.tif', format='tif')
            plt.close()
            return
        else:
            plt.show()    
            return

    


############################################################################
# Functions for PSAM
# 
# 
# 

# Supporting functions for "create_PSAM"

# create list of sequences with Hamming distance 1 at position i, using the alphabet defined in alphabet
def get_changed(sub, i, alphabet='AGCT'):
    """This function has two inputs: sub = sequence, i = position which needs to be altered. Basically it
    creates all sequences which are exactly 1 hamming distance away from sub at position i. 
    E.g get_changed('ACGTTGA', 1) returns ['AAGTTGA', 'AGGTTGA', 'ATGTTGA'].
    
    Important to notice that you can change all_c to for instance single letter AAs if needed. However it will
    only accept single integers as i.
    
    Create list of sequences with Hamming distance 1 at position i, using the alphabet defined in alphabet."""
    
    all_c=set(alphabet)
    other = lambda x : list(all_c.difference(x)) # Lambda function (syntax lambda x: do something with x), other returns list with the other Nucleotides 
    # e.g.other('G') -> ['A','T','C']
    
    return [sub[0:i]+c+sub[i+1:] for c in other(sub[i])]


def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def where_difference(str1, str2):
    """Input two strings with hamming distance 1, return position where difference occurs."""
    return [i for i in range(len(str1)) if str1[i] != str2[i]]


def mutate_hamm_dist1(kmer, alphabet='ACGT'):
    """Similar to the function get_changed, this function creates all possible perturabetion with hamm dist 1
    (dependend on all_c) over all positions in the string.
    
    E.g. mutate_hamm_dist1('ACGT') = ['GCGT','CCGT','TCGT','AAGT','AGGT','ATGT','ACAT',
    'ACCT','ACTT','ACGA','ACGG','ACGC']."""
    
    l = []
    for i in range(len(kmer)):
        l.append(get_changed(kmer,i, alphabet=alphabet))
    l = list(itertools.chain.from_iterable(l))
    return l

# this is obsolete, just use string.replace('CG', 'MG')
def change_CpG_into_MpG(string):
    """This function changes CGs into MGs (of course if the library was methylated)."""
    CG = 'CG'
    change = list(string)
    count = 0
    # search for CG over the entire string (finditer = iterable)
    for match in re.finditer(CG, string):
        count += 1
        # get the position of the C in CG (therefore always index[0]) and change it to M
        change[match.span()[0]] = 'M'
        
    change = "".join(change)
    return change

def powerset(iterable):
    "powerset([1, 2, 3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return [el for el in list(itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1)
    )) if el != ()]

##########################################
# actual PSAM
# 
# Similar to an PWM, a PSAM is a matrix which represent binding affinities 
def create_PSAM(df, column, analysis_mode, fltr=True, return_mutations = False):
    """This function creates a position specific affinity matrix using a tidy dataframe and a modification input.
    Also, pretty neat, methylated CpGs will be displayed as MG to distinguish from nonmethylated CGs. Value = count
    or normalised_counts."""

    if column == 'normalised_counts':
        print('Filtering out infinite values! This might impact the logo if compared to a logo generated with raw counts!')
        df = df[df['normalised_counts'] != np.inf]


    
    
    if analysis_mode == 'pvalue':
        max_kmer = df.nsmallest(1, columns=column)['kmer'].values[0] # extract string
    elif analysis_mode == 'counts':
        max_kmer = df.nlargest(1, columns=column)['kmer'].values[0] # extract string
    else:
        print('Please define analysis mode!') 
        return
    
    


    
    # search for all mutations  with hamming distance 1 with the canonical alphabet
    mutations = mutate_hamm_dist1(max_kmer, alphabet='ACGT')
    
    

    # if at least 1 CG in max kmer, substitute CG with MG and allow hamming distance of greater than 1 to generate psam

    if (max_kmer.count('CG') > 0) | (max_kmer.count('MG') > 0):
        
        if max_kmer.count('CG') > 0:
        # generate kmer with all MGs
            kmer_MG = max_kmer.replace('CG', 'MG')
            mutations.append(kmer_MG)
        elif max_kmer.count('MG') > 0:
        # or a kmer with all CGs
            kmer_CG = max_kmer.replace('MG', 'CG')
            mutations.append(kmer_CG)

        # generate PSAM
        psam = pd.DataFrame(
            {'A':[],
            'C':[],
            'T':[],
            'G':[],
            'M':[]}
        )


        # set psam for max abundant kmer to 1

        for index, letter in enumerate(max_kmer):
                psam.loc[index, letter] = 1

        if analysis_mode == 'pvalue':
            max_kmer_count = df.nsmallest(1, columns=column)['count'].values[0] # extract string
        elif analysis_mode == 'counts':
            max_kmer_count = df.nlargest(1, columns=column)[column].values[0] # extract string
            
    
        for seq in mutations:
            try:
                count = df[df['kmer'] == seq][column].values[0]
            except IndexError:
                count = 0
            
            position_point_mutation = where_difference(max_kmer, seq) # corresponds to the rows in the psam as a list e.g. [1, 4]
            
            for row in position_point_mutation:
                # seq[row] corresponds to 'A, C, G or T' aka the column of the psam
                psam.loc[row, seq[row]] = round(count/max_kmer_count,3)


    #extract all positions where a G appears and allow the prior position to be a M -> to simulate MG context --> M T mimickry, sometimes if there is a TG it actually binds to MG.        
    elif max_kmer.count('G') > 0:
        
        g_starts = [g_match.start() for g_match in re.finditer('G', max_kmer)]
            
        
        #generate all possible subsets (see powerset function)
        all_positions_to_mutate = powerset(g_starts)    
        #iterate over all positions to mutate and change CGs to MGs in all possibilities
        for positions in all_positions_to_mutate:
            
            if positions[0] == 0:
                continue
            
            if len(positions) > 1:
                slices = []
                pointer = 0
                str_to_mutate = max_kmer
                for i, cut_pos in enumerate(positions):
                    slices.append(str_to_mutate[pointer:cut_pos -1] + 'M')
                    pointer = cut_pos 
                    if i == len(positions) - 1:
                        slices.append(str_to_mutate[pointer:])
                mutated_str = ''.join(slices)
            else:
                str_to_mutate = list(max_kmer)
                str_to_mutate[positions[0]-1] = 'M'
                mutated_str = "".join(str_to_mutate)
                
            mutations.append(mutated_str)
        
                 
        # generate PSAM
        psam = pd.DataFrame(
            {'A':[],
            'C':[],
            'T':[],
            'G':[],
            'M':[]}
        )


        # set psam for max abundant kmer to 1

        for index, letter in enumerate(max_kmer):
                psam.loc[index, letter] = 1

    
        if analysis_mode == 'pvalue':
            max_kmer_count = df.nsmallest(1, columns=column)['count'].values[0] # extract string
        elif analysis_mode == 'counts':
            max_kmer_count = df.nlargest(1, columns=column)[column].values[0] # extract string
    
        for seq in mutations:
            try:
                count = df[df['kmer'] == seq][column].values[0]
            except IndexError:
                count = 0
            
            position_point_mutation = where_difference(max_kmer, seq) # corresponds to the rows in the psam as a list e.g. [1, 4]
            
            for row in position_point_mutation:
                # seq[row] corresponds to 'A, C, G or T' aka the column of the psam
                psam.loc[row, seq[row]] = round(count/max_kmer_count,3)
 
    
    else:

        # generate PSAM
        psam = pd.DataFrame(
            {'A':[],
            'C':[],
            'T':[],
            'G':[]
            }
        )

        # set psam for max abundant kmer to 1

        for index, letter in enumerate(max_kmer):
                psam.loc[index, letter] = 1

    
        if analysis_mode == 'pvalue':
            max_kmer_count = df.nsmallest(1, columns=column)['count'].values[0] # extract string
        elif analysis_mode == 'counts':
            max_kmer_count = df.nlargest(1, columns=column)[column].values[0] # extract string

        for seq in mutations:
            try:
                count = df[df['kmer'] == seq][column].values[0]
            except IndexError:
                count = 0
            
            position_point_mutation = where_difference(max_kmer, seq) # corresponds to the rows in the psam as a list e.g. [1, 4]
            
            for row in position_point_mutation:
                # seq[row] corresponds to 'A, C, G or T' aka the column of the psam
                psam.loc[row, seq[row]] = round(count/max_kmer_count,3)


    mutations.append(max_kmer)

    if not return_mutations:
        if not fltr:
            return psam
        
        else:

            for sequence in mutations:
                df = df[df['kmer'] != sequence]
                rt_seq = rev_complement(sequence)
                df = df[df['kmer'] != rt_seq]
            
            return psam, df
    
    elif not fltr:
        return psam, mutations
    
    else:

        for sequence in mutations:
                df = df[df['kmer'] != sequence]
                rt_seq = rev_complement(sequence)
                df = df[df['kmer'] != rt_seq]
        
        return psam, df, mutations








def create_dGG(psam):
    """Takes a psam and yields a dGG energy matrix. for this I quote Foat et al 2006 'Statistical mechanical modeling of genome-wide transcription
        # factor occupancy data by MatrixREDUCE'
    
        #'For each position in the PSAM, the average DDG is calculated. Then, the difference
        #between each individual DDG and the average DDG at that position is
        #computed; the absolute value of this difference is the height of the character
        #representing that nucleotide. If the difference is positive (more favorable
        #than average), the letter is placed above a horizontal black line through the
        #center of the logo. If the difference is negative (less favorable than average)
        #the letter is placed below the black line. Larger letters are stacked on smaller
        #letters moving outward from the black line. The height of the letter can be
        #interpreted as free energy difference from the average in units of RT.'"""
    
    psam = psam + 0.00001
    dGG = np.log(psam)

    if dGG.isnull().values.any():
        # I assume if there are NaN it comes because of ACGTM
        # define rows with and without NaN
        is_NaN = dGG.isnull()
        row_has_NaN = is_NaN.any(axis=1)
        row_no_NaN = dGG[~row_has_NaN]
        rows_with_NaN = dGG[row_has_NaN]

        # split the two dfs and calculate the averages, either with 4 columns or 5 columns
        avg_4 = rows_with_NaN.sum(axis=1) / 4 
        dGG_4 = rows_with_NaN.sub(avg_4, axis = 'rows')
        avg_5 = row_no_NaN.sum(axis=1) / 5
        dGG_5 = row_no_NaN.sub(avg_5, axis = 'rows')

        # merge everything together and concat using the index
        dGG = pd.concat([dGG_4, dGG_5])
        dGG = dGG.sort_index()
        dGG = dGG.fillna(0)
        
        return dGG
    
    else:
        avg = dGG.sum(axis=1) / 4
        #subtract average from entries
        dGG = dGG.sub(avg, axis = 'rows')
        
        return dGG        




###########################################
# now PWM


def estimate_frequencies(df):
    """In order to create an weblogo based on information content, one needs to know the frequencies how often each letter appears(both background and non background).
    This function takes a df in the format created via the function kmer_counting and creates (and returns) a position frequency matrix.
    First, it duplicates rows according to the amount of counts (how often a given kmer is present). Then, it counts how often certain letters appear in certain positions."""

    df = df.loc[df.index.repeat(df['count'])]

    A = pd.DataFrame([list(sequence) for sequence in df['kmer']])
    pfm = A.apply(pd.value_counts)

    pfm = pfm.fillna(0)
    return pfm

def col_wise_propability(column):
    """Supporting function for pfm_to_ppm: Taking one column as input, this function calculates for each entry the propability over the sum of the column.
    Returns a list which can be used new column in position propability matrix."""
    total = sum(column)
    values = []
    for each in column:
        values.append(each/total)
    
    return values

def pfm_to_ppm(df):
    """This function applies on each column the function col_wise_propability and converts therefore a position frequency matrix into a position propability matrix."""
    return df.apply(lambda x: col_wise_propability(x))

def create_PWM(ppm, background = 0.25, bits=2):
    """This function takes as input a ppm and converts it into a PWM. By default, background is set as 0.25, but can be as well adjusted by inputting 
    a different model (in form of a df).
    First, add pseudocounts, then create the log2 of the ratios of df over background. Background can be either float or pd.DataFrame."""

    ppm = ppm + 0.00001
    if bits != 2:
        if type(background) == float:
            pwm = np.log(ppm / background)/ np.log(bits)
    
        elif isinstance(background, pd.DataFrame):
            background = background + 0.00001
            pwm = np.log(ppm / background)/ np.log(bits)

        else:
            print('background must be either a float or a pd.DataFrame')
            return

    else:
        if type(background) == float:
            pwm = np.log2(ppm / background)
    
        elif isinstance(background, pd.DataFrame):
            background = background + 0.00001
            pwm = np.log2(ppm / background)

        else:
            print('background must be either a float or a pd.DataFrame')
            return
    
    return pwm

def estimate_total_IC(ppm, pwm):
    """Calculates total information content using a position propability matrix and position weight matrix. 
    Rsequence(l) := sum(Pij * log ( Pij / Pbg) = ppm * pwm (Schneider et al. 1986)."""

   

    h = ppm * pwm 

    Rseq = h.apply(lambda x: sum(x))

    return Rseq

def convert_PPM_to_IM(ppm, pwm):
    """Converts a position weight matrix into a matrix containing information content for each letter. Rsequence(l) := sum(Pij * log ( Pij / Pbg) = ppm * pwm (Schneider et al. 1986).
    Lettersize is defined as Lettersize = Rsequence * ppm (Schneider and Stephens et al. 1990)."""
    
    

    h = ppm * pwm
    
    Rseq = h.apply(lambda x: sum(x))


    im = Rseq * ppm

    return im

def select_seed(df, column, analysis_mode):
    """Select seed, aka most abundant kmer from dataset and all mutations with hamming distance 1 (or more in case of CG) 
    to create Information logo."""

    if column == 'normalised_counts':
        print('Filtering out infinite values! This might impact the logo if compared to a logo generated with raw counts!')
        df = df[df['normalised_counts'] != np.inf]
    if analysis_mode == 'pvalue':
        max_kmer = df.nsmallest(1, columns=column)['kmer'].values[0] 
    elif analysis_mode == 'counts':
        max_kmer = df.nlargest(1, columns=column)['kmer'].values[0] # extract string
    else:
        print('Please define analysis mode!')
        return
    
    # search for all mutations  with hamming distance 1
    mutations = mutate_hamm_dist1(max_kmer)

    mutations.append(max_kmer)
    
    # if more than 1 CG allow hamming distance of greater than 1 to generate psam

    if max_kmer.count('CG') > 1:
        kmer_MG = max_kmer.replace('CG', 'MG')

        mutations.append(kmer_MG)

    elif max_kmer.count('MG') > 1:
        kmer_CG = max_kmer.replace('MG', 'CG')

        mutations.append(kmer_CG)



    #extract all positions where a G appears and allow the prior position to be a M -> to simulate MG context --> M T mimickry, sometimes if there is a TG it actually binds to MG.        
    elif max_kmer.count('G') > 0:
        
        g_starts = [g_match.start() for g_match in re.finditer('G', max_kmer)]
            
        
        #generate all possible subsets (see powerset function)
        all_positions_to_mutate = powerset(g_starts)    
        #iterate over all positions to mutate and change CGs to MGs in all possibilities
        for positions in all_positions_to_mutate:
            
            if positions[0] == 0:
                continue
            
            if len(positions) > 1:
                slices = []
                pointer = 0
                str_to_mutate = max_kmer
                for i, cut_pos in enumerate(positions):
                    slices.append(str_to_mutate[pointer:cut_pos -1] + 'M')
                    pointer = cut_pos 
                    if i == len(positions) - 1:
                        slices.append(str_to_mutate[pointer:])
                mutated_str = ''.join(slices)
            else:
                str_to_mutate = list(max_kmer)
                str_to_mutate[positions[0]-1] = 'M'
                mutated_str = "".join(str_to_mutate)
                
            mutations.append(mutated_str)

    

    df = df[df['kmer'].isin(mutations)]

    

    return df