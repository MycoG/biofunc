#!/usr/bin/env python3

import argparse
from biofunc.dist import mindist
from biofunc.bed import load_bed, save_bed
import subprocess
from pathlib import Path
from datetime import datetime
import sys
import os

program_name = "abmindist"
description = "calculate minimum distance between two BED files"

def abmindist(A_bed:str, B_bed:str, out_name:str, out_dir:str, dist:int=5_000) -> Path:
    """
    Script to generate minimum distance between two BED files
    
    :param A_bed: Path to BED file
    :type A_bed: str
    :param B_bed: Path to BED file
    :type B_bed: str
    :param out_name: Name for output file (without extension)
    :type out_name: str
    :param out_dir: Output directory
    :type out_dir: str
    :param dist: Maximum distance (in bp) to be counted from a feature in A_bed to B_bed
    :type dist: int
    :return: Absolute Path to mindist file
    :rtype: Path
    """
    #check output directory is valid, try to create if missing
    out_dir_path = Path(out_dir)
    if not out_dir_path.is_dir():
        out_dir_path.mkdir(parents=True, exist_ok=True)

    #create log file
    logpath = out_dir_path/f"{out_name}.{program_name.split('.')[0]}.log"
    open(logpath, "w").close()
    def log_write(str):
        with open(logpath, 'a') as f:
            f.write(f"[{datetime.now()}] {str} \n")
    if __name__ == "__main__":
        log_write(f"COMMAND: {' '.join(sys.argv)}")
        log_write(f"CWD: {os.getcwd()}")

    #check that A_bed, and B_bed are valid BED files, with associated .header
    A_path = Path(A_bed)
    A_header_path =  Path(A_bed+".header")
    B_path = Path(B_bed)
    B_header_path = Path(B_bed+".header")
    for path in [A_path, A_header_path, B_path, B_header_path]:
        if not path.is_file():
            raise Exception(f"path {str(path)} is not a valid path!")
    
    log_write(f"INFO: a_bed -  {A_path.resolve()}")
    with open(A_path, 'r') as f:
        num_lines = sum(1 for _ in f)
        log_write(f"INFO: nrows - {num_lines}")
    log_write(f"INFO: b_bed - {B_path.resolve()}")
    with open(B_path, 'r') as f:
        num_lines = sum(1 for _ in f)
        log_write(f"INFO: nrows -  {num_lines}")

    #declare output files
    outfile = out_dir_path/(out_name+".mindist.bed")
    outfiletmp = out_dir_path/(out_name+".tmp")
    outheadertmp = out_dir_path/(out_name+".tmp.header")

    #perform bedtools window
    command = f"bedtools window -a {str(A_path)} -b {str(B_path)} -w {dist} > {str(outfiletmp)}"
    log_write(f"COMMAND: {command}")
    window_process = subprocess.run(command, shell=True, capture_output=True)
    if window_process.returncode != 0:
        log_write(f"bedtools exited with non-zero status {window_process.returncode}")
        log_write(str(window_process.stdout.decode()))
        log_write(str(window_process.stderr.decode()))
        window_process.check_returncode()

    log_write(f"INFO: {outfiletmp}")
    with open(outfiletmp, 'r') as f:
        num_lines = sum(1 for _ in f)
        log_write(f"INFO: nrows -  {num_lines}")

    #write a headers, b headers, mindist
    log_write(f"Writing tmp headers...")
    with open(outheadertmp, 'w') as f:
        with open(A_header_path, 'r') as g:
            ahead = g.read().strip().split("\n")
            aheaders = ["a_"+x for x in ahead[:3]] + ahead[3:]
        with open(B_header_path, 'r') as g:
            bhead = g.read().strip().split("\n")
            bheaders = ["b_" +x for x in bhead[:3]] + bhead[3:]
        f.write("\n".join(aheaders + bheaders)+"\n")
    log_write(f"Headers written to {str(outheadertmp.resolve())} ! ")

    #calculate mindist
    df = load_bed(str(outfiletmp))
    log_write(f"Calculating mindist...")
    df['mindist'] = df.apply(lambda x: mindist(x['a_start'], x['a_end'], x['b_start'], x['b_end']), axis=1)
    save_bed(df, str(outfile))
    log_write(f"Mindist BED written to {outfile.resolve()}")

    #delete tmp files
    outfiletmp.unlink()
    outheadertmp.unlink()
    log_write(f"Mindist BED tempfiles deleted...")
    log_write(f"Operation Complete! Now go get 'em!")

    return outfile.resolve()

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=program_name,
        description=description
    )
    parser.add_argument("bedfile_A")
    parser.add_argument("bedfile_B")
    parser.add_argument("outname")
    parser.add_argument("outdir")
    parser.add_argument("-d", "--window_distance", default=5000, type=int)
    return parser.parse_args()


def main():
    args:argparse.Namespace = parse_args()
    abmindist(
        A_bed=args.bedfile_A, 
        B_bed=args.bedfile_B,
        out_name=args.outname,
        out_dir=args.outdir,
        dist=args.window_distance)

if __name__ == "__main__":
   main()