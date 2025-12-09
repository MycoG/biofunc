
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.axes import Axes
import pandas as pd

def plot_karyogram(
        ref: pd.DataFrame,
        df: pd.DataFrame,
        ax: Axes = None, 
        height: int = 0.75,
        width: int = None,
        cmap = plt.cm.viridis) -> None:
    """
    plots a karyogram using matplotlib

    :param ref: pandas Dataframe of tabix FASTA reference genome
    :type ref: pd.DataFrame
    :param df: pandas Dataframe of bed file
    :type df: pd.DataFrame
    :param height: height of each of the tracks
    :type height: int
    :param width: number of bp to expand BED file by
    :type width: int
    """

    #TODO: do all the variable refactoring

    fig,ax = plt.subplots(figsize=(16,8), dpi=200)
    height = 0.75
    snp_width = 700_000

    chrom_order = [f"chr{x}" for x in range(1,23)]
    for idx, chr in enumerate(chrom_order):
        df = sup_data[sup_data['chr'] == chr]
        ref = hs1[hs1['NAME'] == chr].iloc[0]

        yloc = idx-height/2

        #plot bars
        ax.broken_barh(
            [[0, ref.LENGTH]],
            (yloc, height),
            color="#FFFFFF", 
            zorder=2
            )

        
        low_psweep_df = df[df['calibrated_smoothed P(sweep)'] < 0.5]
        high_psweep_df = df[df['calibrated_smoothed P(sweep)'] >= 0.5]

        # plot SNPs that are less than 0.5 p(sweep)
        starts = low_psweep_df['start'].to_list()
        intervals = [[x, snp_width] for x in starts]
        ax.broken_barh(
            intervals,
            (yloc, height),
            color=cmap(low_psweep_df['calibrated_smoothed P(sweep)']),
            zorder=2,
            alpha=1
            )
        
        # plot SNPs that are more than 0.5 p(sweep)
        starts = high_psweep_df['start'].to_list()
        intervals = [[x, snp_width] for x in starts]
        ax.broken_barh(
            intervals,
            (yloc, height),
            color=cmap(high_psweep_df['calibrated_smoothed P(sweep)']),
            zorder=3,
            alpha=1
            )

        #plot bars
        ax.broken_barh(
            [[0, ref.LENGTH]],
            (yloc, height),
            color='#00000000', 
            edgecolor='black',
            linestyle="solid",
            linewidth=1.0,
            zorder=4
            )

    ax.set_yticks(range(len(chrom_order)), labels=chrom_order)
    ax.invert_yaxis()

    #grid
    ax.grid(True)

    #hide top, right, bottom spintes
    ax.spines[['top','right']].set_visible(False)
    ax.set_title(f"hs1 karyogram with SwifrSNP annotation\nSNPs expanded by x{snp_width} for visibility")

    #colorbar
    sm = ScalarMappable(cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="P(sweep)", aspect=50)

    plt.show()

def plot_chrom():
    """
    plot a chromosome using a BED file
    """
    pass