from .het import _calc_hets

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