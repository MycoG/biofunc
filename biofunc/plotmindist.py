#!/usr/bin/env python3

import argparse
from biofunc.bed import load_bed
from pathlib import Path
from datetime import datetime
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

program_name = "plotmindist"
description = "plot features according to minimum distance column, and log some summary stats"

def _plot(df:pd.DataFrame, 
          x:str, 
          y:str, 
          ran, 
          outfile, 
          suptitle:str=None, 
          axtitle:str=None,
          window_size:int=50,
          ) -> None:
    
    fig, ax = plt.subplots(dpi=500)

    #scatterplot
    ax.scatter(df['mindist'], df[y], alpha=0.2)

    sorted_df = df.sort_values("mindist")
    
    #plotting by slicing mindist seems to show similar results to plotting above
    def roll_func(x):
        slice:pd.DataFrame = sorted_df[y][ (sorted_df['mindist'] >= x) & (sorted_df['mindist'] < x+window_size) ]
        return slice.mean()
    sorted_df['rolled'] = sorted_df['mindist'].apply(roll_func)
    ax.plot(sorted_df['mindist'], sorted_df['rolled'], label=f"new - window={window_size}bp", color="green")

    #plot rolling mean by mindist
    rolling_df = sorted_df[["mindist",y]].rolling(window_size, center=True, on="mindist").mean()
    ax.plot(rolling_df["mindist"], rolling_df[y], label=f"original - window={window_size}", color="red")
    

    #config
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_xlim(-ran, ran)
    ax.set_ylim(0, 1)
    ax.set_title(axtitle, fontsize=8)
    plt.legend(loc='upper right')
    fig.suptitle(suptitle, fontsize=10)

    #save figure
    fig.savefig(outfile)

    return

def _gen_quant(df:pd.DataFrame, col:str):
    q_val = 0.01
    high_quant = df[col].quantile(1-q_val)
    low_quant = df[col].quantile(q_val)
    return q_val, low_quant, high_quant

def _get_quant_df(df:pd.DataFrame, col:str, low_quant:float, high_quant:float):
    high_filt_df = df[df[col] >= high_quant]
    low_filt_df = df[df[col] <= low_quant]
    mid_filt_df = df[ (df[col] < high_quant) & (df[col] > low_quant) ] \
        .sample(int(np.average([len(high_filt_df), len(low_filt_df)])))
    return low_filt_df, mid_filt_df, high_filt_df

def mindist_plot(input_bed:str, cols:str, ranges:str, out_name:str, out_dir:str, pltlabel:str, rolling_window_size:int=50, quant_col:str="calibrated_smoothed P(sweep)"):
    
    #create output directory if not valid
    out_dir:Path = Path(out_dir)
    out_dir_plots = out_dir / "plots"
    if not out_dir.is_dir():
        out_dir.mkdir(parents=True, exist_ok=True)
        out_dir_plots.mkdir(parents=True, exist_ok=True)

    #setup logs
    log = open(out_dir/f"{program_name}.{out_name}.log", "w")
    def log_write(str):
        log.write(f"[{datetime.now()}] {str} \n")
    if __name__ == "__main__":
        log_write(f"COMMAND: {" ".join(sys.argv)}")

    #check if input bed is valid
    input_bed:Path = Path(input_bed)
    if not input_bed.is_file():
        log_write(f"ERR: File {input_bed} is not valid!")
        raise Exception(f"File {input_bed} is not valid!")
    
    #load input bed
    bedfile:pd.DataFrame = load_bed(str(input_bed))
    cols:list = _get_cols(cols, list(bedfile.columns), log_write)
    ranges:list = _get_ranges(ranges, log_write)

    #check pltlabel
    if pltlabel == None:
        pltlabel = ""

    #generate quantiles
    q_val, low_qval, high_qval = _gen_quant(bedfile, quant_col)
    low_df, mid_df, high_df = _get_quant_df(bedfile, quant_col, low_qval, high_qval)
    log_write(f"INFO: q_col - {quant_col}")
    log_write(f"INFO: q >= {high_qval} - n={len(high_df)}")
    log_write(f"INFO: {low_qval} >= q >= {high_qval} - n={len(mid_df)} (sampled)")
    log_write(f"INFO: q <= {low_qval} - n={len(low_df)}")

    quant_tuples = [
        ("low-0.01", f"{quant_col}<={low_qval:.4f}", low_df),
        ("mid-0.98", f"{low_qval:.3f}<{quant_col}<{high_qval:.4f}", mid_df),
        ("top-0.01", f"{quant_col}>={high_qval:4f}", high_df),
        ("all", f"all qval", bedfile.sample(len(mid_df)))
    ]
    
    #plot by each quantile, column, and range
    for label, qval_str, df in quant_tuples:
        for col in cols:
            for ra in ranges:

                out_dir_range = out_dir_plots / (str(ra)+"kb")
                out_dir_range.mkdir(parents=True, exist_ok=True)
                outfile = out_dir_range / f"{label}_{col}_{ra}.png"

                #plot
                # suptitle = f"{col} | {label} | {ra}kb"
                suptitle = " | ".join([col, pltlabel, (str(ra)+"kb from central SNP")])
                axtitle = f"{label} | {qval_str}"
                _plot(df=df, x="mindist", y=col, ran=ra, suptitle=suptitle, axtitle=axtitle, outfile=outfile, window_size=rolling_window_size)

                print(f"wrote {label}, col: {col}, range: {ra}, path: {str(outfile)}")
                
                # within set range
                # log number of number of rows to A 
                # and number of rows to B
                range_df:pd.DataFrame = df[ (df['mindist'] >= - ra) & (df['mindist'] <= ra)]
                range_groupby = range_df.groupby(by=['a_chrom', 'a_start'])

                log_write(f"INFO: wrote {label}, col: {col}, range: {ra}, path: {str(outfile)}")
                log_write(f"INFO: num of features in A - {range_groupby.ngroups}")
                log_write(f"INFO: num features in B - {len(range_df)}")

    log_write("Operation Complete! Now go get 'em! ")
    return

def _get_cols(cols:str, df_cols: list, log_write) -> list:
    "get columns from program input, else ask for columns interactively"
    if cols == None:
        #ask interactively
        while True:
            print("please input columns as a comma separated list!")
            print(f"valid columns are {df_cols} :")
            cols:list = input().replace(", ", ",").replace(" ,", ",").split(",")
            col_diff = set(cols).difference(set(df_cols))
            if col_diff == set():
                break
    else:
        # cols is set
        cols:list = cols.replace(", ", ",").replace(" ,", ",").split(",")
        if len(cols) == 0:
            log_write(f"ERR: Number of cols is 0")
            log_write(f"ERR: col list - {cols}")
            raise Exception(f"Number of columns is 0")
    log_write(f"INFO: cols selected - {cols}")
    print(f"cols selected - {cols}")
    return cols

def _get_ranges(ranges:str, log_write) -> list:
    if ranges == None:
        #ask interactively
        while True:
            print("please input distance from features (ex. 5000 is 5kb away from central feature)")
            ranges:list = input().replace(", ", ",").replace(" ,", ",").split(",")

            #check if all numbers are digits
            if all([x.isdigit() for x in ranges]):
                ranges = [int(x) for x in ranges]
                break
    else:
        ranges:list = ranges.replace(", ", ",").replace(" ,", ",").split(",")
        if len(ranges) == 0:
            log_write(f"ERR: Number of ranges is 0")
            raise Exception(f"Number of ranges is 0")
        if all([x.isdigit() for x in ranges]):
            ranges = [int(x) for x in ranges]
        else:
            raise Exception(f"Range contains non-digit")
    log_write(f"INFO: ranges selected - {ranges}")
    print(f"ranges selected - {ranges}")
    return ranges

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=program_name,
        description=description
    )
    parser.add_argument("input_BED")
    parser.add_argument("out_name")
    parser.add_argument("out_dir")

    #optional args - if not set the program asks interactively
    parser.add_argument("-c", "--cols", help="comma-separated list of columns to plot")
    parser.add_argument("-r", "--ranges", help="comma-separated list of ranges to plot")
    parser.add_argument("-l", "--label")
    parser.add_argument("-w", "--rolling_window", help="rolling window size (in bp)", default=50)
    parser.add_argument("-q", "--quantile_column", default="calibrated_smoothed P(sweep)")
    return parser.parse_args()

def main():
    args:argparse.Namespace = parse_args()
    
    window_size = int(args.rolling_window)
    mindist_plot(
        input_bed=args.input_BED,
        cols=args.cols,
        ranges=args.ranges,
        out_name=args.out_name,
        out_dir=args.out_dir,
        pltlabel=args.label,
        rolling_window_size=window_size,
        quant_col=args.quantile_column
    )

if __name__ == "__main__":
    main()