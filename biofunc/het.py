import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def convert_het_bed(path, out=None) -> None:
    """
    converts a '.het' file generated using calc_het into a BED file.

    :param path: The path to a *.het file
    :type path: str
    :param out: Output file path
    :type out: str
    """
    het = pd.read_csv(path, sep='\t')
    het['END'] = het['POS'] + 1
    het_filt = het[['CHROM', 'POS', 'END', 'O(Het)', 'E(Het)' ]]

    #
    if out != None:
        outpath = out
    else:
        outpath = path+".bed"

    #write BED file
    het_filt.to_csv(outpath, sep="\t", header=False, index=False)
    print(f"bed file written to {outpath}.bed")

    #write HEADER file
    with open(outpath+".header", 'w') as f:
        f.write("\n".join(het_filt.columns)) 
        print(f"header written to {outpath}.bed.header")

    return

def calc_het(path) -> None:
    """
    calculates heterozygosity using 

    :param path: A path to a VCF file.
    :type path: str
    """
    #TODO: add output paramater to control filename, location output

    def __calc_hets(line:str):
        """
        Calculate observed and expected heterozygosity per TR

        @param line: A non-header line from a VCF file
        """

        num_alt = len(line[4].split(",")) #line[4].count(",") + 1
        samples = line[9:]
        n_samples = len(samples)

        #create a dictionary to keep count of all the alleles
        allele_dict = { str(n):0 for n in range(num_alt + 1)} #number of alt + 1 reference

        num_het = 0
        for sample in samples:
            sample = sample.split(":")

            # get genotype field
            gt = sample[0].replace("|", "/").split("/")
            
            #
            if len(gt) > 2: 
                continue

            if gt[0] != gt[1]:
                num_het += 1

            #add to allele dict
            allele_dict[ str(gt[0]) ] += 1
            allele_dict[ str(gt[1]) ] += 1

        est_hom = 0 
        for k,fi in allele_dict.items():                         #summation of alleles
            est_hom += ( fi )  ** 2                              # pi = (count of allele i) / (2 * number of indiv)
        
        het_obs = num_het / len(samples)                        # obs_heterozygosity = (# heterozygous individuals / N individuals )
        het_est = 1 - (est_hom / ( 2 * n_samples ) ** 2 )        # est_heterozygosity = (1 - est homozygosity)

        return het_obs, het_est
    
    f = open(path, 'r')
    g = open(path+".het", "w")

    for line in f:
        n_header_cols = 8

        #skip header lines
        if line.startswith("##"):
            continue
        #get header line
        elif line.startswith("#"):
            header = line.strip()[1:].split("\t")
            g.write("\t".join(header[:n_header_cols]) + "\tO(Het)\tE(Het)\n")

        else:
            line = line.strip().split("\t")
            het_obs, het_est = __calc_hets(line)

            g.write("\t".join([ str(t) for t in line[:n_header_cols] ] + [str(het_obs), str(het_est)]) + "\n")

    f.close()
    g.close()
    
    return