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

import logomaker

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



######################################################################################
# For scripting and calculations

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


def rev_complement_seq(seq):
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

def where_CG(seq):
    return [i for i in range(len(seq) - 1) if seq[i:i+2] == 'CG']

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
        l.append(single_point_mut_seq(kmer,i, alphabet=alphabet))
    l = list(itertools.chain.from_iterable(l))
    return l


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

def calculate_ratios(kmer_df):
    """This function takes kmer df and 'normalizes' (i.e. creates ratios eluted/input) 
    the eluted kmers with the input. """

    kmer_e = kmer_df[kmer_df['status'] == 'eluted']
    # write all kmers as a regex pattern
    pattern = '|'.join(kmer_e['kmer'])
    # Use .str.contains with regex=True to filter kmer_df for all eluted kmers
    filtered_kmer = kmer_df[kmer_df['kmer'].str.contains(pattern, regex=True)]

    # iterate over the filtered kmer dataframe and divide the eluted count by input 
    # counts (methylated and unmethylated separately)

    m_dict = {}
    nm_dict = {}
    for k, df in filtered_kmer.groupby('kmer'):
        # extract count of methylated kmer in input
        im_count = df[(df['mod'] == 'methl') & (df['status'] == 'input')]['count'].values
        # generate ratio and save as dictionary
        m_dict[k] = df[df['mod'] == 'methl']['count'].values / im_count
        # extract count of nonmethylated kmer in input
        inm_count = df[(df['mod'] == 'nonmethl') & (df['status'] == 'input')]['count'].values
        # generate ratio
        nm_dict[k] = df[df['mod'] == 'nonmethl']['count'].values / inm_count

    # store everything in a new df
    new_df = pd.DataFrame.from_dict(m_dict, orient='index', columns=['input_methl', 'eluted_methl']).T.replace([], np.nan)
    new_df = new_df.T.reset_index()

    new_df2 = pd.DataFrame.from_dict(nm_dict, orient='index', columns=['input_nonmethl', 'eluted_nonmethl']).T.replace([], np.nan)
    new_df2 = new_df2.T.reset_index()

    #add CG information
    final_df = pd.merge(new_df, new_df2, on='index')
    final_df['CpG'] = ['CG' in final_df['index'][i] for i in range(len(final_df))]

    return final_df



#######################################################################################################
# For all psam/ppm matrix related stuff

def rev_complement(df, extended_alphabet=True):
    """Takes a df where columns are named 'ACGT', or 'ACGTmg' and inverts the matrix."""
    if extended_alphabet:
        #rev complement
        rev_colnames = ['T','G','C','A', 'g', 'm']
        df_rev = df.iloc[::-1]
        df_rev.columns = rev_colnames
        df_rev = df_rev.reset_index(drop=True)
        df_rev = df_rev[['A','C','G','T','m','g']]

        return df_rev
    
    else:
        rev_colnames = ['T','G','C','A']
        df_rev = df.iloc[::-1]
        df_rev.columns = rev_colnames
        df_rev = df_rev.reset_index(drop=True)
        df_rev = df_rev[['A','C','G','T']]

        return df_rev

def calculate_stats2(dGG):
    """Calculate the letter sizes of mg and CG form dGG logos, but take only values into account if they are > 0
    (that is if mg or CG are above 0, otherwise it indeicates that both CG and mg are disadvantages for the TF-DNA interaction).
    Calcualte the difference between the pairs, and calculate how large the letter size is in relation to the overall information
    stored at the given position in the logo (expressed in percent). Returns a stats_df comprising average CG, and mg values, 
    the delta between those two, the ratio CG/mg, and the percentages of C, G, m and g.
    """
    CG = []
    mg = []
    ratio_CGmg = []
    delta_s = []
    per_c = []
    per_G = []
    per_m = []
    per_g = []
    
    limit = len(dGG)
    euclidean_norms = dGG.apply(lambda row: row.abs().sum(), axis=1) # euclidean norm over the psam positions
    top_val = dGG.values.max()
    for i in range(limit-1):
        C = dGG.loc[i, 'C']
        G = dGG.loc[i+1, 'G']
        m = dGG.loc[i, 'm']
        g = dGG.loc[i+1, 'g']
        
        delta = (m - C) + (g - G) # calculate delta as suggested by reviewer, if positive, mg favoured, if negative CG.
        delta_s.append(delta)

        if (C > 0) & (G > 0):
            t = (C + G) / 2
            CG.append(t)
        else: 
            CG.append('NA')

        if(m > 0) & (g > 0):
            tm = (m + g) / 2
            mg.append(tm)
        else:
            mg.append('NA')

        try:
            ratio = CG[i] / mg[i]
        except TypeError:
            if (type(CG[i]) == str) & (type(mg[i]) == str):
                ratio = 'not applicable'
            elif type(CG[i]) == str:
                ratio = 'mg available'
            elif type(mg[i]) == str:
                ratio = 'CG available'
        ratio_CGmg.append(ratio)
        
        
        # calculating percentages not based on top value but on the euclidean norm (the absolute spread of the psam)
        # -1 shows directionality, negative value is not beneficial for binding
        try:
            percent_c = C/euclidean_norms[i]
        except TypeError:
            percent_c = 'not applicable'
        per_c.append(percent_c)

        try:
            percent_G = G/euclidean_norms[i+1]
        except TypeError:
            percent_G = 'not applicable'
        per_G.append(percent_G)
        
        try:
            percent_m = m/euclidean_norms[i]
        except TypeError:
            percent_m = 'not applicable'
        per_m.append(percent_m)
        
        try:
            percent_g = g/euclidean_norms[i+1]
        except TypeError:
            percent_g = 'not applicable'
        per_g.append(percent_g)
   

    stats_df = pd.DataFrame({'CG': CG, 
                             'mg': mg, 
                             'ratio CG/mg': ratio_CGmg, 
                             'C%topvalue': per_c,
                             'G%topvalue': per_G,
                             'm%topvalue': per_m,
                             'g%topvalue': per_g,
                             'delta_mg-CG': delta_s
                            })
    return stats_df

def single_point_mut_seq(sub, i, alphabet='AGCT'):
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

#Do I need powerset?
def powerset(iterable):
    "powerset([1, 2, 3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return [el for el in list(itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1)
    )) if el != ()]

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

#Do I need estimate_freqs?
def estimate_frequencies(df):
    """In order to create an weblogo based on information content, one needs to know the frequencies how often each letter appears(both background and non background).
    This function takes a df in the format created via the function kmer_counting and creates (and returns) a position frequency matrix.
    First, it duplicates rows according to the amount of counts (how often a given kmer is present). Then, it counts how often certain letters appear in certain positions."""

    df = df.loc[df.index.repeat(df['count'])]

    A = pd.DataFrame([list(sequence) for sequence in df['kmer']])
    pfm = A.apply(pd.value_counts)

    pfm = pfm.fillna(0)
    return pfm

def _col_wise_propability(column):
    """Supporting function for pfm_to_ppm: Taking one column as input, this function calculates for each entry the propability over the sum of the column.
    Returns a list which can be used new column in position propability matrix."""
    total = sum(column)
    values = []
    for each in column:
        values.append(each/total)
    
    return values

def ddG_to_psam(df):
    """Converts a ddG (Probound output) into psam."""
    mean_values = df.max(axis=1)



    df_m = df.add(-mean_values, axis='rows')



    # ProBound uses ln
    psam = np.e**(df_m)
    return psam

def pfm_to_ppm(df):
    """This function applies on each column the function col_wise_propability and converts therefore a position frequency matrix into a position propability matrix."""
    df = df.T
    return df.apply(lambda x: _col_wise_propability(x)).T

# DOn't think I need create_PWM
def create_PWM(ppm, background = 0.25, bits=2):
    """This function takes as input a ppm and converts it into a PWM. By default, background is set as 0.25, but can be as well adjusted by inputting 
    a different model (in form of a df).
    First, add pseudocounts, then create the log2 of the ratios of df over background. Background can be either float or pd.DataFrame."""

    ppm = ppm + 0.00000001
    if bits != 2:
        if type(background) == float:
            pwm = np.log(ppm / background)/ np.log(bits)
    
        elif isinstance(background, pd.DataFrame):
            background = background + 0.00000001
            pwm = np.log(ppm / background)/ np.log(bits)

        else:
            print('background must be either a float or a pd.DataFrame')
            return

    else:
        if type(background) == float:
            pwm = np.log2(ppm / background)
    
        elif isinstance(background, pd.DataFrame):
            background = background + 0.00000001
            pwm = np.log2(ppm / background)

        else:
            print('background must be either a float or a pd.DataFrame')
            return
    
    return pwm

# Don't think i need this neither
def estimate_total_IC(ppm, pwm):
    """Calculates total information content using a position propability matrix and position weight matrix. 
    Rsequence(l) := sum(Pij * log ( Pij / Pbg) = ppm * pwm (Schneider et al. 1986)."""

   

    h = ppm * pwm 

    Rseq = h.apply(lambda x: sum(x))

    return Rseq

#same
def convert_PPM_to_IM(ppm, pwm):
    """Converts a position weight matrix into a matrix containing information content for each letter. Rsequence(l) := sum(Pij * log ( Pij / Pbg) = ppm * pwm (Schneider et al. 1986).
    Lettersize is defined as Lettersize = Rsequence * ppm (Schneider and Stephens et al. 1990)."""
    
    

    h = ppm * pwm
    
    Rseq = h.apply(lambda x: sum(x))


    im = Rseq * ppm

    return im

def ddg_to_consensus(ddg_df: pd.DataFrame) -> str:
    """
    Generate a consensus sequence from a -ΔΔG matrix.

    Parameters:
    - ddg_df: pandas DataFrame, shape (L, 4), where rows are positions and
              columns are nucleotides ('A', 'C', 'G', 'T').
              Higher values represent more favorable binding (e.g., -ΔΔG).

    Returns:
    - consensus: str, the consensus sequence.
    """
    consensus = ''
    for _, row in ddg_df.iterrows():
        consensus += row.idxmax()
    return consensus

def find_stabilizing_CG(ddg_df):
    """
    Identify positions where a CG dinucleotide is stabilizing for TF-DNA binding.
    A CG dinucleotide is considered stabilizing if:
    - The ΔΔG for C at position i > 0
    - The ΔΔG for G at position i+1 > 0

    Parameters
    ----------
    ddg_df : pd.DataFrame
        A DataFrame with columns ['A', 'C', 'G', 'T'] and rows as positions.

    Returns
    -------
    List[int]
        List of positions i where a CG dinucleotide (C at i, G at i+1) is stabilizing.
    """
    positions = []

    for i in range(len(ddg_df) - 1):  # stop at second-last row
        c_val = ddg_df.loc[i, 'C']
        g_val = ddg_df.loc[i + 1, 'G']
        if c_val > 0 and g_val > 0:
            positions.append(i)

    return positions


###################################################################################################################################
# For CG statistics

def estimate_nucleotide_distribution(
    df, 
    seq_col, 
    methylation_col=None, 
    methylation_BC=None, 
    apply_methylation_encoding=False
):
    """
    Estimates nucleotide distribution from sequences in a DataFrame, optionally selecting only one library 
    (methylated or unmethylated) and applying CG → m/g substitution for methylated sequences.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing DNA sequences.
    seq_col : str
        Column name with 24N DNA sequences.
    methylation_col : str or None, optional
        Column name containing molecular barcodes identifying methylation state.
    methylation_BC : str or None, optional
        The barcode value identifying the library species to analyze (e.g. 'AGTA' or 'GAAT').
    apply_methylation_encoding : bool, optional
        If True and methylation_BC is selected, CG dinucleotides will be converted to 'm'/'g'.
        Default is True.

    Returns
    -------
    pd.Series
        Normalized nucleotide frequencies (A, C, G, T, m, g as applicable).
    """
    all_nucs = []

    # Filter by methylation barcode if specified
    if methylation_col and methylation_BC:
        df = df[df[methylation_col] == methylation_BC]

    for seq in df[seq_col].dropna():
        seq = seq.upper()

        if apply_methylation_encoding and methylation_col and methylation_BC:
            # Replace CG with m/g for methylated sequences
            i = 0
            while i < len(seq) - 1:
                if seq[i:i+2] == 'CG':
                    all_nucs.extend(['m', 'g'])
                    i += 2
                else:
                    all_nucs.append(seq[i])
                    i += 1
            if i == len(seq) - 1:
                all_nucs.append(seq[-1])
        else:
            all_nucs.extend(seq)

    counts = collections.Counter(all_nucs)
    total = sum(counts.values())

    return pd.Series({k: v / total for k, v in sorted(counts.items())})

def estimate_dinucleotide_distribution(
    df, 
    seq_col, 
    methylation_col=None, 
    methylation_BC=None, 
    apply_methylation_encoding=False
    ):
    """
    Estimates dinucleotide distribution (including methylation-aware 'mg').

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing DNA sequences.
    seq_col : str
        Name of the column with DNA sequences.
    methylation_col : str or None, optional
        Column indicating molecular barcode for methylation (optional).
    methylation_BC : str or None, optional
        Barcode value identifying methylated sequences.

    Returns
    -------
    pd.Series
        Normalized dinucleotide frequencies (e.g., AA, AT, mg, gA, etc.)
    """
    dinucs = []

    for _, row in df.iterrows():
        seq = row[seq_col].upper()

        mod_seq = []
        i = 0
        if apply_methylation_encoding and methylation_col and methylation_BC and row.get(methylation_col) == methylation_BC:
            while i < len(seq) - 1:
                if seq[i:i+2] == 'CG':
                    mod_seq.extend(['m', 'g'])
                    i += 2
                else:
                    mod_seq.append(seq[i])
                    i += 1
            if i == len(seq) - 1:
                mod_seq.append(seq[-1])
        else:
            mod_seq = list(seq)

        for i in range(len(mod_seq) - 1):
            dinucs.append(mod_seq[i] + mod_seq[i + 1])

    counts = collections.Counter(dinucs)
    total = sum(counts.values())
    return pd.Series({k: v / total for k, v in counts.items()})

def compute_kl_divergence(ppm: pd.DataFrame, background: pd.Series) -> pd.Series:
    """
    Computes KL divergence at each position in a Position Probability Matrix (PPM),
    relative to a given background distribution. Supports both standard (A, C, G, T)
    and extended (e.g., m, g) alphabets.

    Parameters
    ----------
    ppm : pd.DataFrame
        Position Probability Matrix (rows: positions, columns: nucleotides).
        Columns can include A, C, G, T, and optionally m, g.

    background : pd.Series
        Background nucleotide distribution. Index must match all columns of `ppm`.
        For example: ['A', 'C', 'G', 'T', 'm', 'g'] if using an extended alphabet.

    Returns
    -------
    pd.Series
        KL divergence per position (indexed like `ppm`).
    """
    # Ensure background includes all bases in the PPM
    bases = [base for base in ppm.columns if base in background.index]
    ppm = ppm[bases].copy()
    background = background[bases]

    kl_div = []

    for i, row in ppm.iterrows():
        kl = 0
        for base in bases:
            p = row[base]
            q = background[base]

            if p > 0 and q > 0:
                kl += p * np.log2(p / q)
        kl_div.append(kl)

    return pd.Series(kl_div, index=ppm.index, name='KL_divergence')

def jensen_shannon_divergence(P: pd.DataFrame, Q: pd.DataFrame) -> pd.Series:
    """Per-position Jensen-Shannon divergence, using two pd.Dataframes of the same length!"""
    P = P.replace(0, 1e-10)
    Q = Q.replace(0, 1e-10)
    M = 0.5 * (P + Q)
    return 0.5 * (P * np.log2(P / M)).sum(axis=1) + 0.5 * (Q * np.log2(Q / M)).sum(axis=1)

def crop_consensus(con_series, threshold=0.1):
    """Crops a pd.Series based on threshold.

    Parameters
    ----------
    con_series : pd.Series
    threshold : float, optional

    Returns
    -------
    pd.Series
    """
    mask = con_series > threshold
    if not mask.any():
        print('None of the values surpass the threshold, returning original series.')
        return con_series  

    first_idx = mask.idxmax()
    last_idx = mask[::-1].idxmax()
    return con_series.loc[first_idx:last_idx]

def resolve_cg_js_div(CG, js_div):
    """
    Resolves divergent positions between two sequences while accounting for CG dinucleotides.

    This function ensures that CG dinucleotides are treated as linked:
    - If neither the C nor G in a CG pair is in js_div, the C is kept.
    - If only the C is in js_div, the C is kept.
    - If only the G is in js_div, the G is kept (instead of the C).
    - If both are in js_div, both positions are kept.
    
    Additionally, all divergent positions that are not part of a CG dinucleotide are kept.

    Parameters:
    ----------
    CG : set of int
        Positions of the first base ('C') in each CG dinucleotide (0-based index).
    js_div : iterable of int
        Positions where two sequences differ, based on Jensen-Shannon divergence.

    Returns:
    -------
    List[int]
        Sorted list of final positions to keep, incorporating both CG logic and other divergent positions.
    """
    CG = set(CG)  # ensure it's a set
    js_div = set(js_div)  # ensure it's a set

    final = set()

    for c_pos in CG:
        g_pos = c_pos + 1
        c_in = c_pos in js_div
        g_in = g_pos in js_div

        if c_in and not g_in:
            final.add(c_pos)
        elif not c_in and g_in:
            final.add(g_pos)
        elif c_in and g_in:
            final.update([c_pos, g_pos])
        else:
            final.add(c_pos)  # neither in js_div, but still keep C

    # Add js_div positions that are not part of a CG dinucleotide (i.e., not Gs that follow a C)
    cg_g_positions = {c + 1 for c in CG}
    remaining_js_div = js_div - cg_g_positions
    final.update(remaining_js_div)

    return sorted(final)

def optimal_consensus_alignment(cropped_m, cropped_um, dGGm, dGGum, max_shift=2):
    """
    Try extending or trimming both methylated and unmethylated ΔΔG windows
    to align consensus sequences with the minimal Hamming distance.

    Parameters:
    - cropped_m, cropped_um: DataFrames with .index from KL-cropping
    - dGGm, dGGum: full ΔΔG DataFrames
    - max_shift: int, how far (±) to shift/extend/trim from original range

    Returns:
    - best_idx_m, best_idx_um: pd.Index ranges for dGGm and dGGum
    """

    range_m = list(range(cropped_m.index[0], cropped_m.index[-1] + 1))
    range_um = list(range(cropped_um.index[0], cropped_um.index[-1] + 1))

    best_m, best_um = range_m, range_um
    min_hamming = float('inf')

    # try combinations of shifts: (start_shift, end_shift)
    for shift_m_start in range(-max_shift, max_shift + 1):
        for shift_m_end in range(-max_shift, max_shift + 1):
            m_start = range_m[0] + shift_m_start
            m_end = range_m[-1] + shift_m_end
            idx_m = list(range(m_start, m_end + 1))

            if m_start < dGGm.index[0] or m_end > dGGm.index[-1] or m_start >= m_end:
                continue

            for shift_um_start in range(-max_shift, max_shift + 1):
                for shift_um_end in range(-max_shift, max_shift + 1):
                    um_start = range_um[0] + shift_um_start
                    um_end = range_um[-1] + shift_um_end
                    idx_um = list(range(um_start, um_end + 1))

                    if um_start < dGGum.index[0] or um_end > dGGum.index[-1] or um_start >= um_end:
                        continue

                    # generate consensus sequences for both ranges
                    cons_m = ddg_to_consensus(dGGm.loc[idx_m, :])
                    cons_um = ddg_to_consensus(dGGum.loc[idx_um, :])

                    # they must be equal length to compute Hamming
                    if len(cons_m) != len(cons_um):
                        continue

                    dist = hamming_distance(cons_m, cons_um)

                    if dist < min_hamming:
                        min_hamming = dist
                        best_m, best_um = idx_m, idx_um

    return pd.Index(best_m), pd.Index(best_um)

####################################################################################################################################
# For all plots

def plot_messages(messages, ax):
    """
    Display messages in a textbox-like format on the given axis.
    """
    ax.axis('off')
    message_text = "\n".join(messages)
    ax.text(0, 1, message_text, fontsize=9, verticalalignment='top', family='monospace')

def get_stacked_ylim(*dfs):
    """
    Given one or more ddG matrices, return the y-axis limits 
    needed to display all letter stacks properly.
    """
    max_stacks = []
    min_stacks = []
    for df in dfs:
        df = df.fillna(0)  # in case there are NaNs
        max_stacks.append(df[df > 0].sum(axis=1).max())
        min_stacks.append(df[df < 0].sum(axis=1).min())
    return min(min_stacks), max(max_stacks)

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

def make_logo(logomatrix, ax=None, title=None, ylim=None, ylabel="$-\Delta \Delta G$ (kcal/mol)"):
    """
    Plots a logo on the given Axes. Most commonly used with ddG/psam logos. Suboptimal for ppms due to dimensions.

    Parameters:
    - logomatrix: pd.DataFrame compatible with logomaker.Logo
    - ax: Optional matplotlib Axes object. If None, a new figure is created.
    - ylim: Optional tuple (ymin, ymax) to set y-axis limits.
    - ylabel: Label the y-axis.
    
    Returns:
    - ax: The matplotlib Axes containing the logo.
    """
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=[10, 6])
        created_fig = True
    
    logo = logomaker.Logo(logomatrix,
                        shade_below=0.5,
                        ax=ax,
                        fade_below=0.5,
                        color_scheme={'A':'#66a61e', 'C':'#7570b3','G':'#ffc809','T':'#d95f02','m':'#a6cee3'}
                        )
    # style using Logo methods
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left'], visible=True)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    logo.ax.set_ylabel(ylabel, labelpad=-1)
    logo.ax.xaxis.set_ticks_position('none')
    logo.ax.xaxis.set_tick_params(pad=-1)

    if ylim is not None:
        logo.ax.set_ylim(ylim)

    if created_fig:
        plt.tight_layout()
        plt.show()
    
    if title:
        logo.ax.set_title(title, fontsize=10)
    return ax  

def ppm_logo(df, ax=None, title=None, ylim=None, ylabel="Probability"):
    """Like make_logo, but better for ppm matrices."""
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=[7, 2])
        created_fig = True
    logo = logomaker.Logo(df,
                          ax=ax,
                          shade_below=0.5,
                          fade_below=0.5,
                          color_scheme={'A':'#66a61e', 'C':'#7570b3','G':'#ffc809','T':'#d95f02','m':'#a6cee3'},
                          stack_order = 'small_on_top' 
                          )

    # style using Logo methods
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0)
        

        # style using Axes methods
    logo.ax.set_ylabel(ylabel, labelpad=-1)
    logo.ax.xaxis.set_ticks_position('none')
    logo.ax.xaxis.set_tick_params(pad=-1)

    if ylim is not None:
        logo.ax.set_ylim(ylim)
    
    if title:
        logo.ax.set_title(title, fontsize=10)
    
    if created_fig:
        plt.tight_layout()
        plt.show()
        return
    
    return ax    

def plot_delta_mg_CG_context_aware(stats_df, threshold=0.1, figsize=(9, 4), ax=None):
    """
    Plot Δ(mg − CG) values across positions in a dinucleotide logo matrix,
    showing only values that meet two criteria:

    1. The position must have valid CG or mg context (i.e., not 'NA' in 'CG' or 'mg' columns).
    2. The dinucleotide basepairs CG or mg must have a relative height ≥ `threshold`
       compared to the position's total signal (i.e., Euclidean norm).

    Invalid positions are grayed out with 'x' markers.

    Parameters
    ----------
    stats_df : pandas.DataFrame
        A DataFrame created by `calculate_stats2`, containing Δ(mg − CG) values and the context
        and relative height of base letters across dinucleotide positions.

    threshold : float, optional
        The minimum relative height (e.g. 0.1 = 10%) to consider a position meaningful. Default is 0.1.

    figsize : tuple, optional
        Size of the matplotlib figure (used only if ax is None). Default is (9, 4).

    ax : matplotlib.axes.Axes or None
        Optional. If provided, the plot will be drawn on this axis. Otherwise, a new figure is created.

    Returns
    -------
    matplotlib.axes.Axes
        The axis object used for plotting.
    """

    x = stats_df.index
    y = stats_df['delta_mg-CG']

    # Determine context availability
    has_context = (stats_df['CG'] != 'NA') | (stats_df['mg'] != 'NA')

    # Compute average letter height for CG and mg
    cg_contrib = stats_df[['C%topvalue', 'G%topvalue']].mean(axis=1)
    mg_contrib = stats_df[['m%topvalue', 'g%topvalue']].mean(axis=1)
    has_height = (cg_contrib > threshold) | (mg_contrib > threshold)

    valid = has_context & has_height

    # Setup figure and axis if not provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True

    # Plotting
    ax.plot(x[valid], y[valid], marker='o', linestyle='-', color='black', label='Δ(mg − CG)')
    ax.scatter(x[~valid], y[~valid], marker='x', color='gray', label='Filtered out')

    ax.axhline(0, color='lightgray', linestyle='--')

    ax.set_title("Δ(mg − CG) per Dinucleotide Position", fontsize=12)
    ax.set_xlabel("Letter Position (i)", fontsize=10)
    ax.set_ylabel("Δ(mg − CG)", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis='y', linestyle=':', linewidth=0.5)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax

def kmer_scatterplots(ratio_df, coords=None, ax=None, figsize=(5.5, 5.5)):
    """
    Plot methylated vs. unmethylated DNA elution ratios for enriched kmers,
    as used in Gralak et al.

    Parameters
    ----------
    ratio_df : pandas.DataFrame
        DataFrame containing eluted_methl, eluted_nonmethl, CpG, and fill_plot columns.
    
    coords : dict, optional
        Dictionary of regression line coordinates for 'with_CG' and/or 'no_CG'.
        Each key should map to {'x': ..., 'y': ...}. Default is None.

    ax : matplotlib.axes.Axes, optional
        Matplotlib axis to plot on. If None, a new figure will be created.

    figsize : tuple, optional
        Used only if ax is None. Default is (5.5, 5.5).

    Returns
    -------
    matplotlib.axes.Axes
        The axis object used for plotting.
    """
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

    ratios = ratio_df  # alias

    x = ratios['eluted_methl']
    y = ratios['eluted_nonmethl']

    # Prepare colors
    edgecolors = ratios['CpG'].map({True: darkred, False: lightblack})
    facecolors = [
        edgecolors.iloc[i] if sig else 'None'
        for i, sig in enumerate(ratios['fill_plot'])
    ]

    # Create figure/axis if needed
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True

    # Scatter plot
    ax.scatter(x=x, y=y, facecolors=facecolors, edgecolors=edgecolors, rasterized=True)

    # Regression lines
    if coords:
        if 'with_CG' in coords:
            ax.plot(coords['with_CG']['x'], coords['with_CG']['y'], color=darkred)
        if 'no_CG' in coords:
            ax.plot(coords['no_CG']['x'], coords['no_CG']['y'], color=lightblack)

    # Aesthetics
    ax.grid(visible=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_aspect('equal', adjustable='box')


    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())


    lower_limit = min_val - (max_val * 0.05)
    upper_limit = max_val * 1.05

    ax.set_xlim(lower_limit, upper_limit)
    ax.set_ylim(lower_limit, upper_limit)

    ax.set_xlabel('methylated DNA', fontfamily='sans-serif', fontsize=10, fontstyle='italic')
    ax.set_ylabel('unmethylated DNA', fontfamily='sans-serif', fontsize=10, fontstyle='italic')

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax

def plot_kl_divergence(kl_series, threshold_minor=0.1, threshold_major=0.3, ax=None, title="KL divergence per position"):
    """
    Plot KL divergence values per position with reference threshold lines.

    Parameters
    ----------
    kl_series : pd.Series
        KL divergence per position (index = position, values = KL values).

    threshold_minor : float, optional
        Lower threshold to visualize minor information content. Default is 0.1.

    threshold_major : float, optional
        Upper threshold to visualize high information content. Default is 0.3.

    ax : matplotlib.axes.Axes, optional
        If provided, the plot is drawn on this axis. Otherwise, creates a new figure.

    title : str, optional
        Title of the plot. Default is "KL divergence per position".

    Returns
    -------
    None
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3))

    # Plot KL divergence as line
    ax.plot(kl_series.index, kl_series.values, color='black', linewidth=1.2, marker='o')

    # Reference thresholds
    ax.axhline(threshold_minor, linestyle=':', color='gray', label=f'{threshold_minor} (minor info)')
    ax.axhline(threshold_major, linestyle='--', color='red', label=f'{threshold_major} (major info)')

    # Axis labels and title
    ax.set_xlabel("Position", fontsize=10)
    ax.set_ylabel("KL Divergence", fontsize=10)
    ax.set_title(title, fontsize=11)

    # Format ticks and spines
    ax.tick_params(axis='both', labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Show legend if standalone
    if ax is None:
        ax.legend(frameon=False, fontsize=8)

    if ax is None:
        plt.tight_layout()
        plt.show()
