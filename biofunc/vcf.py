from .het import _calc_hets
import gzip

head_idx = {
    "CHROM": 0,
    "POS": 1,
    "ID": 2,
    "REF": 3,
    "ALT": 4,
    "QUAL": 5,
    "FILTER": 6,
    "INFO": 7,
    }

def _get_info_dict(info_col) -> dict:
        info_dict = {}
        info = info_col.split(";")
        for item in info:
            item:str
            split_item = item.split("=")
            info_dict[split_item[0]] = split_item[1]
        return info_dict

def vcf_to_bed(input:str, output:str, header:bool=False, expand_info:bool=False) -> None:
    """
    Convert BED file to VCF

    :param input: path to a VCF file
    :type input: str, path
    :param output: path to a BED file
    :type output: str, path
    :param header: Write header into bed file?
    :type header: bool
    :param expand_info: Expands INFO column into separate columns instead of just one
    :type expand_info: bool
    """

    input_file = open(input, "r")
    output_file = open(output, "w")
    headers_written = False

    #write mandatory BED headers
    if header or expand_info:
        output_file.write("chrom\tchromStart\tchromEnd")

    for line in input_file:
        line = line.strip()

        if line.startswith("#"):
            continue

        else:
            line = line.split("\t")
            info:dict = _get_info_dict(line[head_idx["INFO"]])
            chrom_end = info.pop("END")

            if expand_info:
                if not headers_written:
                    info_headers = info.keys()
                    output_file.write("\t" + "\t".join(info_headers) + "\n")
                    headers_written = True

            het = _calc_hets(line=line)
            het = [str(x) for x in het]
            # Write chrom, chromStart, chromEnd, and columns if enabled
            out_line = "\t".join(line[:2] + [chrom_end] + [v for k,v in info.items() if expand_info] + het) + "\n"
            output_file.write(out_line)
        

    input_file.close()
    output_file.close()

    return

def _calculate_het(line):
    line = line.split("\t")
    gt_lines = line[9:]

    num_alt_alleles = len(line[4].split(","))
    allele_dict = { str(k):0 for k in range(0,num_alt_alleles+1) }

    num_het = 0
    num_indiv = 0
    for i, genotype in enumerate(gt_lines):
        genotype:str

        alleles = genotype.split("|")
        #skip more /less than 2 alleles
        if len(alleles) != 2:
            continue
        
        a1, a2 = alleles

        #check if contains missing genotypes
        #if it does, do not include in calculation
        if ((a1 == ".") or (a2 == ".")):
            #skip if either allele is missing
            continue
        else:
            num_indiv += 1
            if a1 != a2:
                num_het += 1

            allele_dict[a1] += 1
            allele_dict[a2] += 1
    
    if num_indiv == 0:
        return (None, None)

    Ho = num_het / num_indiv
    Homo_e = 0
    for allele, freq in allele_dict.items():
        Homo_e += ( freq / (2*num_indiv) ) ** 2

    He = 1-Homo_e
    return Ho, He
    
def calc_het(path, out):
    f = open(path, 'r')
    g = open(out, "w")

    for line in f.readlines():
        line = line.strip()
        
        if line.startswith("#"):
            continue
        else:
            hets = _calculate_het(line)
            chrom, pos = line.split("\t")[:2]
            g.write("\t".join([str(x) for x in [chrom, pos, pos]] + [str(x) for x in hets])  + "\n")
            
    
    f.close()
    g.close()

    #write header
    with open(out+".header", 'w') as f:
        f.write("\n".join(["chrom", "start", "end", "O(Het)", "E(Het)"]))

def calc_het_gz(path, out):
    #get num lines in path
    print("finding num of lines...")
    with gzip.open(path, 'rb') as f:
        num_lines = sum(1 for _ in f)
    print(f"num of lines: {num_lines}")

    open(out, "w").write("")
    
    with gzip.open(out, 'at') as g:
        with gzip.open(path, 'rt') as f:
            print("opened file...")
            for i, line in enumerate(f):
                print(f"{i} / {num_lines} processed.", end="\r")
                line = line.strip()
                if line.startswith("#"):
                    continue
                else:
                    hets = _calculate_het(line)

                    #if hets is invalid, return empty strings
                    if hets == (None, None):
                        hets = ["", ""]

                    chrom, pos = line.split("\t")[:2]
                    g.write("\t".join([str(x) for x in [chrom, pos, pos]] + [str(x) for x in hets])  + "\n")
            
    #write header
    with open(out+".header", 'w') as f:
        f.write("\n".join(["chrom", "start", "end", "O(Het)", "E(Het)"]))