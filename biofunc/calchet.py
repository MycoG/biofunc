#!/usr/bin/env python3

import argparse
from pathlib import Path
import gzip
from datetime import datetime
import os
import sys
import time
import numpy as np

program_name = "calchet.py"
description = "script to calculate observed and expected heterozygosity"

def calc_het(input_vcf, output_bed, output_dir, compressed=False):
    input_vcf:Path = Path(input_vcf)
    output_dir:Path = Path(output_dir)
    output_bed:Path = output_dir / f"{output_bed}.bed"
    output_header:Path = Path(str(output_bed)+".header")
    log_path:Path = output_dir / f"{program_name.split('.')[0]}.log"

    log_write = _validate_input(input_vcf, output_dir, output_bed, log_path)

    if compressed:
        input_file = gzip.open(input_vcf, 'rt')
    else:
        input_file = open(input_vcf, 'r')

    output_file = open(output_bed, 'w')

    log_write("looping over input file...")
    for line in input_file:
        line = line.strip()

        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line[1:].split("\t")
            with open(output_header, "w") as f:
                f.write("chrom\nstart\nend\nO(Het)\nE(Het)\n")
            log_write(f"INFO: header written to: {output_header}")
        #calculate het
        else:
            line=line.split("\t")
            het_obs, het_est = _calc_het_line(line)
            output_file.write("\t".join([ str(t) for t in line[:2] ] + [str(line[1])] + [str(het_obs), str(het_est)]) + "\n")

    input_file.close()
    output_file.close()
    log_write("Operation complete! Now go get 'em!")
    pass

def _calc_het_line(line:str):
    """Calculate observed and expected heterozygosity per TR"""

    num_alt = len(line[4].split(",")) #line[4].count(",") + 1
    samples = line[9:]
    # n_samples = len(samples)
    valid_samples = 0

    #create a dictionary to keep count of all the alleles
    allele_dict = { str(n):0 for n in range(num_alt + 1)} #number of alt + 1 reference

    num_het = 0
    for sample in samples:
        sample = sample.split(":")

        # genotype field per sample
        # format example: [0, 1]
        gt = sample[0].replace("|", "/").split("/")
        
        if (gt[0] == ".") or (gt[1] == "."):
            continue

        # Check if genotype(gt) is heterozygous
        if gt[0] != gt[1]:
            num_het += 1

        #add to allele dict to get allele count per TR
        allele_dict[ str(gt[0]) ] += 1
        allele_dict[ str(gt[1]) ] += 1
        valid_samples += 1

    if valid_samples == 0:
        return np.nan, np.nan

    # Observed Heterozygosity
    het_obs = num_het / valid_samples                      # obs_heterozygosity = (# heterozygous individuals / N individuals )


    #Estimated Heterozygosity
    est_hom = 0 
    for allele ,freq in allele_dict.items():                #summation of alleles squared
        est_hom += ( freq )  ** 2                           # pi = (count of allele i) / (2 * number of indiv)
    het_est = 1 - (est_hom / ( 2 * valid_samples ) ** 2 )       # est_heterozygosity = (1 - est homozygosity)
    
    return het_obs, het_est

def _validate_input(input_vcf, output_dir, output_bed, log_path):

    # check if output bed directory is valid
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True, exist_ok=True)

    #create log file
    open(log_path, "w").close()
    def log_write(str):
        with open(log_path, 'a') as f:
            f.write(f"[{datetime.now()}] {str} \n")
    if __name__ == "__main__":
        log_write(f"COMMAND: {' '.join(sys.argv)}")
        log_write(f"CWD: {os.getcwd()}")


    # check if vcf is valid path
    if not input_vcf.is_file():
        log_write(f"Input VCF is not valid! path: {input_vcf}")
        raise Exception("Input VCF is not valid")

    return log_write
    
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=program_name,
        description=description
    )
    parser.add_argument("input_vcf")
    parser.add_argument("out_name")
    parser.add_argument("out_dir")
    parser.add_argument("-c", "--compressed", action="store_true")
    return parser.parse_args()

def main():
    args:argparse.Namespace = parse_args()
    start_time = time.time()
    calc_het(
        input_vcf=args.input_vcf,
        output_bed=args.out_name,
        output_dir=args.out_dir,
        compressed=args.compressed
             )
    end_time = time.time()
    print(f"took {end_time-start_time:.2f}s")

if __name__ == "__main__":
    main()