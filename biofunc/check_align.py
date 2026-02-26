
from Bio import SeqIO
from pathlib import Path
from pysam import FastaFile
import multiprocessing as mp
from argparse import ArgumentParser
import pandas as pd


def check_align(input_bed_path:str, ref_fasta_A, ref_fasta_B, outdir):
    bed = pd.read_csv(input_bed_path, sep="\t", names=["a_chr", "a_start", "a_end", "b_chr", "b_start", "b_end"], chunksize=100)
    for df in bed:
        print(df)
        quit()

    pass

def get_args():
    parser = ArgumentParser(
        prog="check_align.py"
    )
    parser.add_argument("input_bed")
    parser.add_argument("fasta_ref_A")
    parser.add_argument("fasta_ref_B")
    parser.add_argument("-d", "--output-dir")
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    check_align(
        input_bed=args.input_bed,
        ref_fasta_A=args.ref_fasta_A,
        ref_fasta_B=args.ref_fasta_B,
        outdir=args.output_dir
    )