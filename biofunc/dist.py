import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from pathlib import Path

def mindist(x_start:int, x_end:int, y_start:int, y_end:int) -> int:
    """
    calculate minimum bp distance between two features X and Y

    :param x_start: Start position of feature X
    :type x_start: int
    :param x_end: End position of feature X
    :type x_end: int
    :param y_start: Start position of feature Y
    :type y_start: int
    :param y_end: End position of feature Y
    :type y_end: int
    :return: The minimum distance between X and Y
    :rtype: int
    """
    a = y_start - x_start
    b = y_end - x_start 
    c = y_start - x_end
    d = y_end - x_end
    min = np.min(np.abs([a,b,c,d]))
    for dist in [a, b, c, d]:
        if np.abs(dist) == min:
            return dist


def plot_mindist(df: pd.DataFrame, col:str, out:str,
                range:int=50_000, 
                suptitle:str=None,
                scale='linear',
                window_size=1000,
                ) -> None:
    """
    Plot distance plots. Requires min-dist column.

    :param df: pandas dataframe with "min-dist" column
    :type df: pd.DataFrame
    :param col: column to plot on y-axis
    :type col: str
    :param out: output filepath
    :type out: str
    :param range: distance from x-axis to bounds
    :type range: int
    :param suptitle: title of the figure
    :type suptitle: str
    :param scale: scale for the y-axis of the plot
    :type scale: 'linear' OR 'log'
    :param window_size: size of sliding windows
    :type window_size: int
    """
    fig, ax = plt.subplots()
    d = df.sort_values('min-dist')

    ax.scatter(d['min-dist'], d[col], alpha=0.1)

    # plot centered rolling mean
    d['rolling-mean'] = d[col].rolling(window_size, center=True).mean()
    ax.plot(d['min-dist'], d['rolling-mean'], color='red', label="rolling-mean", linewidth=2)

    #plot global mean
    global_mean = d[col].mean()
    ax.hlines(global_mean, xmin=-500_000, xmax=500_000, linestyles='dotted', linewidth=2, color='indigo', label=f'mean = {global_mean:.3f}')

    ax.set_ylabel(col)
    ax.set_yscale(scale)
    ax.set_xlabel('dist (bp)')
    ax.set_xlim(-range, range)

    size = len( d[ (d['min-dist'] > -range) & (d['min-dist'] < range) ] )
    ax.set_title(f"{range}bp from SNP | window_size={window_size} | n={size}")
    ax.legend(loc="upper right")
    fig.suptitle(suptitle)
    plt.savefig(out)