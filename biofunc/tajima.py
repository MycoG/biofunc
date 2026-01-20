from allel import iter_vcf_chunks, windowed_tajima_d, GenotypeArray

def windowed_tajima(vcf, window_size, step, outfile, min_sites=3):
    """Iterate over VCF chunks and write Tajima's D as BED"""
    #open up vcf, make sure it is sorted
    fields, samples, headers, it = iter_vcf_chunks(input=vcf, fields=["variants/CHROM", "variants/POS", "calldata/GT"], region="chr1")
    
    #leftover numpy array
    left_chrom = None
    left_pos = None
    left_gt = None
    for chunk in it:
        chunk_dict = chunk[0]
        chromosomes = chunk_dict['variants/CHROM']
        positions = chunk_dict['variants/POS']
        gt_array = chunk_dict['calldata/GT']

        #if there are leftover genotypes, thenadd them to gt_array for next tajima_d
        
        #save end of gt_array into variable when using step sizes
        #length of gt_leftover should be window_size - 1
        left_gt = gt_array[-window_size-1:]

        #allele counts
        ac = GenotypeArray(gt_array).count_alleles()

        #perform windowed tajima
        D, windows, counts = windowed_tajima_d(positions, size=window_size, step=step, )


        print(chunk)
        quit()
    
    #loop over w
    pass