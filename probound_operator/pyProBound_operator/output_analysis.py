import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pyProBound import ProBoundModel


def __get_slides(sequence, size):
    slides = [sequence[i:i + size] for i in range(len(sequence) - size + 1)]
    return slides


def plot_kmers_scores(sequences,
                      eluted_counts,
                      model,
                      labels=None, ax=None,
                      color_label_map=None,
                      size=1,
                      counts_percentile=1-10e-5, alpha=0.75
                      ):
    """
    Plots relationship between the k-mer counts and ProBound scores
    :param sequences: probe sequences
    :param eluted_counts: eluted counts of the probe sequences
    :param model: fitted ProBoundModel
    :param labels: labels of probe sequences for plotting
    :param ax: axis to plot on. If None, a new figure is created
    :param color_label_map: label to matplotlib color map
    :param size: point size
    :param counts_percentile: percentile for top cuttof for counts percentile
    :param alpha: alpha for plotting
    :return: fig, ax OR only ax if predefined
    """
    input_ax = ax
    df = pd.DataFrame(sequences)
    df.columns = ["sequence"]
    df["scores"] = model.score_binding_mode_scores(df["sequence"].values,
                                                   score_format="profile",
                                                   profile_aggregate="max")
    df["scores"] = df["scores"].str[0]  # only one binding mode is used
    psam_size = len(df["sequence"].iloc[0]) - len(df["scores"].iloc[0]) + 1
    df["slides"] = df["sequence"].apply(lambda sq: __get_slides(sq, psam_size))
    df["eluted"] = eluted_counts
    if labels is not None:
        df["label"] = labels
        groupers = ["slides", "label", "scores"]
    else:
        groupers = ["slides", "scores"]
    data = df.explode(["slides", "scores"])
    gdf = data.groupby(by=groupers).agg({"eluted": np.sum}).reset_index()

    if input_ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        ax = input_ax

    if labels is None:
        ax.scatter(x=gdf["scores"], y=gdf["eluted"], s=size, alpha=alpha)
    elif color_label_map is None:
        ax.scatter(x=gdf["scores"], y=gdf["eluted"], c=gdf["label"], s=size, alpha=alpha)
    else:
        for key, val in color_label_map.items():
            mask = gdf["label"] == key
            ax.scatter(x=gdf.loc[mask, "scores"], y=gdf.loc[mask, "eluted"],
                       c=val, s=size, label=key, alpha=alpha)

    eluted_lim = np.quantile(gdf["eluted"], counts_percentile)
    ax.set_ylim(bottom=-0.5, top=eluted_lim + 2.5)

    # plot description
    ax.legend(loc='upper right')
    ax.set_xlabel("ProBound scores")
    ax.set_ylabel(f"Eluted k-mer counts, k={psam_size}")
    ax.set_title("Correlation between\nk-mer count and ProBound score.")

    # return result
    if input_ax is None:
        return fig, ax
    else:
        return ax


def plot_comparison(sequences, counts_x, counts_y, kmersize=8,
                    axes_titles=("counts x", "counts y"), color_label_map=None,
                    labels=None, input_ax=None, size=1.0, alpha=0.75):
    """
    Using a sliding window, each sequence is transformed fo kmers of specified size.
    Number of observations is calculated for each kmer (separately for each label).
    Observations are then plotted.
    :param alpha: alpha for plotting
    :param size: point size
    :param input_ax: axis to plot on
    :param color_label_map: label to color map (dictionary)
    :param ax: axis to plot to. If None (default), a new figure will be made.
    :param sequences: iterable of sequences of arbitrary length.
    :param counts_x: counts to plot on x axis
    :param counts_y: counts to plot on y axis
    :param kmersize: size of kmers
    :param axes_titles: how to set labels of x and y axes
    :param labels: iterable of labels for plotting, same size as sequence
    :return: figure, ax OR only ax if ax is not None
    """
    df = pd.DataFrame(sequences)
    df.columns = ["sequence"]
    df["counts_x"], df["counts_y"] = counts_x, counts_y
    df["slides"] = df["sequence"].apply(lambda sq: __get_slides(sq, kmersize))
    if labels is not None:
        df["label"] = labels
        groupers = ["slides", "label"]
    else:
        groupers = ["slides"]

    data = df.explode(["slides"])
    gdf = data.groupby(by=groupers).agg({
        "counts_x": np.sum, "counts_y": np.sum
    }).reset_index()

    if input_ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        ax = input_ax

    if labels is None:
        ax.scatter(x=gdf["counts_x"], y=gdf["counts_y"], s=size, alpha=alpha)
    elif color_label_map is None:
        ax.scatter(x=gdf["counts_x"], y=gdf["counts_y"], c=gdf["label"], s=size, alpha=alpha)
    else:
        for key, val in color_label_map.items():
            mask = gdf["label"] == key
            ax.scatter(x=gdf.loc[mask, "counts_x"], y=gdf.loc[mask, "counts_y"],
                       c=val, s=size, label=key, alpha=alpha)

    ax.legend(loc='upper right')
    ax.set_xlabel(axes_titles[0])
    ax.set_ylabel(axes_titles[1])
    ax.set_title(f"Correlation between\nk-mer counts for k={kmersize}")

    # return result
    if input_ax is None:
        return fig, ax
    else:
        return ax


def methylated_smileseq_analysis(dataframe,
                                 model,
                                 binding_mode=1,
                                 figsize=(8, 4),
                                 mg_coding="mg",
                                 meth_labels=("methylated", "unmethylated")
                                 ):
    """
    Does a complete analysis of mSMiLe-seq output data. Plot a correlation with
    :param meth_labels: labels distinguishing between "methylated" and "unmethylated"
    :param mg_coding: encoding of the methylated CpG
    :param figsize: size of the output figure
    :param model: a path to fit file from ProBound OR a ProBoundModel instance
    :param binding_mode: binding mode id
    :param dataframe: a pandas dataframe with columns:
        "sequence": probe sequence, methylated CGs are translated to <mg_coding>
        "input": input count
        "eluted": eluted count
        "label": label of the dataset: <meth_labels>
    :return: figure, axes
    """

    if type(model) == str:
        model = ProBoundModel(model, fitjson=True)
        model.select_binding_mode(binding_mode)

    meth_label, unmeth_label = meth_labels

    colormaps = {
        f"{unmeth_label}": "dodgerblue",
        f"{meth_label}, CG in probe": "crimson",
        f"{meth_label}, unaffected": "orange"}
    cgprobe_colormap = {"unaffected": "orange", "CG in probe": "teal"}

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=figsize)

    # plot kmers: methylated vs. unmethylated, with labels on whether probe has mg/CG
    dataframe["seq_original"] = dataframe["sequence"].str.replace(mg_coding, "CG", regex=False)
    pivotted = pd.pivot_table(dataframe,
                              values=["eluted"],
                              index=["seq_original"],
                              columns=["label"],
                              aggfunc=np.sum
                              ).fillna(0).reset_index()
    pivotted.columns = ["_".join(x) for x in pivotted.columns]
    pivotted["CG_probe"] = "unaffected"
    pivotted.loc[pivotted["seq_original_"].str.contains("CG"), "CG_probe"] = "CG in probe"
    plot_comparison(pivotted["seq_original_"],
                    pivotted["eluted_methylated"],
                    pivotted["eluted_unmethylated"],
                    labels=pivotted["CG_probe"],
                    axes_titles=("methylated", "unmethylated"),
                    color_label_map=cgprobe_colormap,
                    input_ax=axes[0],
                    size=0.5,
                    kmersize=6
                    )

    # plot kmers vs scores
    methyl_probe_mask = dataframe["sequence"].str.contains(mg_coding)
    methyl_unm = (dataframe["label"] == meth_label) & ~methyl_probe_mask
    methyl_meth = (dataframe["label"] == meth_label) & methyl_probe_mask
    dataframe.loc[methyl_unm, "label"] = dataframe.loc[methyl_unm, "label"] + ", unaffected"
    dataframe.loc[methyl_meth, "label"] = dataframe.loc[methyl_meth, "label"] + ", CG in probe"

    plot_kmers_scores(dataframe["sequence"],
                      dataframe["eluted"],
                      model,
                      labels=dataframe["label"],
                      ax=axes[1],
                      color_label_map=colormaps
                      )
    plt.tight_layout()
    return fig, axes
