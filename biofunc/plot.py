import pandas as pd
from matplotlib.axes import Axes
import matplotlib.pyplot as plt



def plot_mean(
        df, 
        col, 
        ax) -> None:
    """
    plot mean of all data as horizontal line

    """
    #plot global mean
    global_mean = df[col].mean()
    ax.hlines(global_mean, xmin=-500_000, xmax=500_000, linestyles='dotted', linewidth=2, color='indigo', label=f'mean = {global_mean:.3f}')
    return

def plot_rolling_mean(
        df:pd.DataFrame,
        order_by:str,
        on:str,
        ax:Axes,
        size:int=10) -> None:
    """
    plot rolling mean of a dataframe

    :param df: pandas dataframe
    :type df: pd.DataFrame
    :param order_by: column to sort by
    :type order_by: str
    :param on: column to calculate rolling mean
    :type on: str
    :param ax: matplotlib axes
    :type ax: optional, matplotlib.axes.Axes
    :param size: Number of data points to calculate rolling mean by
    :type:
    """
    df = df.sort_values(order_by)
    # plot centered rolling mean
    df['rolling-mean'] = df[on].rolling(size, center=True).mean()
    if ax != None:
        ax.plot(df['min-dist'], df['rolling-mean'], color='red', label=f"rolling-mean of {on}", linewidth=2)
    else:
        plt.plot(df['min-dist'], df['rolling-mean'], color='red', label=f"rolling-mean of {on}", linewidth=2)
    return

def plot_interpolated_mean(
        df: pd.DataFrame,
        x: str,
        y: str,
        ):
    """
    
    :param:
    """
    #groupby some column, and calculate the mean
    data = df.copy().sort_values(group_by)
    groupby = data.groupby(group_by).mean()
    #interpolate the calculated mean
    return

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