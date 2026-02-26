from pysam import VariantFile, VariantRecord
from typing import Tuple
from argparse import ArgumentParser
import numpy as np

program_name = "calchet"
description = "calculate observed and expected heterozygosity per variant in VCF"

def calchet(input_vcf, output_bed):
    vcf = VariantFile(input_vcf, 'r')

    #write file
    with open(output_bed, 'w') as f:

        for rec in vcf.fetch():
            chrom = rec.chrom
            start = rec.start
            stop = rec.stop
            
            #maybe parallelize this?
            het_obs, het_est = _calchet_record(rec)

            format = [chrom, start, stop, het_obs, het_est]
            f.write("\t".join([str(x) for x in format]) + "\n")

    #write header
    with open(str(output_bed)+".header", 'w') as f:
        f.write("chrom\nstart\nstop\nO(Het)\nE(Het)\n")

    return

def _calchet_record(rec:VariantRecord) -> Tuple[float, float]:
    #for each record
    num_alt = len(rec.alts) if rec.alts != None else 0
    n_samples = len(rec.samples)

    #create allele dict
    allele_dict = { str(n):0 for n in range(num_alt + 1) }
    
    #loop over samples
    num_het = 0
    num_indiv = 0
    for sample in rec.samples.values():
        a1, a2 = sample.allele_indices

        #skip if a1 or a2 are None
        if (a1 == None) or (a2 == None):
            continue

        #check if heterozygous
        if a1 != a2 :
            num_het += 1
        
        #add to allele dict
        allele_dict[ str(a1) ] += 1
        allele_dict[ str(a2) ] += 1
        num_indiv += 1

    if num_indiv == 0:
        return (np.nan, np.nan)

    het_obs = num_het / num_indiv

     #Estimated Heterozygosity
    est_hom = 0 
    for allele ,freq in allele_dict.items():                #summation of alleles squared
        est_hom += ( freq )  ** 2                           # pi = (count of allele i) / (2 * number of indiv)
    het_est = 1 - (est_hom / ( 2 * n_samples ) ** 2 )       # est_heterozygosity = (1 - est homozygosity)
    
    return (het_obs, het_est)

def main():
    args = parse_args()
    calchet(
        input_vcf=args.input_vcf,
        output_bed=args.output_bed
        )

def parse_args():
    parser = ArgumentParser(prog=program_name, description=description)
    parser.add_argument("input_vcf")
    parser.add_argument("output_bed")
    return parser.parse_args()

if __name__ == "__main__":
    main()