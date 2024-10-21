import pandas as pd
from pyfastaq.sequences import file_reader as fa_reader


def read_sequences(filename):
    """
    Read sequence file.
    :param filename: path to the input file in fasta, fastq, gzipped fasta or gzipped fastq format.
    If the file is gzipped, the filename suffix is expected to be ".gz".
    :return: pandas dataframe with columns [header, sequence]
    """
    reader = fa_reader(filename)  # reads fasta, fastq, gzipped fasta/q -- includes sanitation
    df = pd.DataFrame([(entry.id, entry.seq) for entry in reader])
    df.columns = ["header", "sequence"]
    return df


def pick_reads(seq_dataframe, nmax):
    """
    Randomly pick nmax reads from a dataframe with sequences.
    :param seq_dataframe: pandas dataframe with columns [header, sequence]
    :param nmax: number of reads to pick
    :return: filtered dataframe of size nmax
    If seq_dataframe is smaller than nmax, seq_dataframe will not be changed (no upsampling).
    """
    if len(seq_dataframe) <= nmax:
        return seq_dataframe
    # choice = np.random.choice(seq_dataframe.index, size=nmax)
    return seq_dataframe.sample(nmax)


def __check_equal_size(*seq_dataframes):
    req_size = len(seq_dataframes[0]["sequence"].values[0])

    for df in seq_dataframes:
        sizes = df["sequence"].str.len().unique()
        if (len(sizes) != 1) | (sizes[0] != req_size):
            raise ValueError("Sequnces of different sizes found: " + str(sizes))


def __calculate_counts(seq_dataframe):
    # should be super efficient even on noninteger data according to the pandas documentation
    counts = seq_dataframe["sequence"].value_counts().reset_index()
    counts.columns = ["sequence", "count"]
    counts["count"] = counts["count"].astype(int)
    return counts


def build_count_table(*seq_dataframes, output_filename=None, gzip=False):
    """
    Build a count table.
    :param gzip: whether to gzip the output count table.
    If True, output_filename will be extended by a ".gz" suffix (if not present).
    :param seq_dataframes: dataframes with columns [header, sequence] to combine in a count table.
    Order of the dataframes is kept. The count columns will have the same order.
    Sequences must be same size.
    :param output_filename: name of the file to save the count table (in valid input format for ProBound).
    If None (default), count table will not be saved.
    :return: count table dataframe
    """
    if len(seq_dataframes) == 0:
        raise ValueError("Nothing to build count table from.")
    # check for equal size of kmers
    __check_equal_size(*seq_dataframes)
    # checking passed -- all sequences are of equal length

    # calculate sequence occurences
    # merge into a count table
    count_table = None
    count_cols = []
    for i, df in enumerate(seq_dataframes):
        current = __calculate_counts(df)
        if count_table is None:
            count_table = current
            count_cols.append("count")
        else:
            count_table = pd.merge(count_table,
                                   current,
                                   on="sequence",
                                   how="outer",
                                   suffixes=("", f"_{i}"),
                                   ).fillna(0)
            count_cols.append(f"count_{i}")

    for count_col in count_cols:
        count_table[count_col] = count_table[count_col].astype(int)

    # order -- just to be sure
    count_cols = sorted(count_cols)
    count_table = count_table[["sequence", *count_cols]]

    # save count table
    if output_filename is not None:
        if gzip:
            if not output_filename.endswith(".gz"):
                output_filename = output_filename + ".gz"
                # infer compression in to_csv func will take care of the rest

        count_table.to_csv(output_filename, sep="\t", header=None, index=False)
    return count_table
