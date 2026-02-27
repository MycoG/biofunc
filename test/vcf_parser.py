import gzip
from pathlib import Path
from typing import Tuple, List

class VCF():
    def __init__(self, path:str, compressed:bool=False):
        self.path:str = path
        self.compressed:bool = compressed
        self.samples: list
        self._len : int = 0

        self._META = {}
        self.PEDIGREE : List[dict] = []
        self.SAMPLE : List[dict] = []
        self.CONTIG : dict = {}
        self.PEDIGREE_DB = None
        self.ALT : List[dict] = []
        self.FORMAT : dict = {}
        self.INFO : dict = {}
        self.FILTER : dict = {}

        input_file = gzip.open(path, 'rt') if self.compressed else open(path, 'r')
        for line in input_file:
            line : str
            if line.startswith("##"):
                line = line.removeprefix("##")
                self._handle_meta(line)
                self._len -= 1
            elif line.startswith("#"):
                line = line.removeprefix("#").strip().split("\t")
                self._handle_header(line)
                self._len -= 1
            else:
                self._len += 1
        input_file.close()

    @property
    def len(self):
        """number of variants within the VCF file"""
        return self._len

    def _handle_header(self, cols:list):
        #first 8 columns should always be the same - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
        #columns 9 is FORMAT, so col 10 and beyond must be sample data
        try:
            self.samples = cols[9:]
        except IndexError:
            self.samples = None

    def _handle_meta(self, line:str):
        line = line.strip()
        # file format
        if line.startswith("file"):
            pass
        # INFO fields
        if line.startswith("INF"):
            line = line.removeprefix("INFO=")
            line_dict = self._handle_xml_fmt(line)
            id = line_dict.pop("ID")
            self.INFO[id] = line_dict
        # Filter fields
        if line.startswith("FIL"):
            line = line.removeprefix("FILTER=")
            line_dict = self._handle_xml_fmt(line)
            id = line_dict.pop("ID")
            self.FILTER[id] = line_dict
        # Individual Format
        if line.startswith("FOR"):
            line = line.removeprefix("FORMAT=")
            line_dict = self._handle_xml_fmt(line)
            id = line_dict.pop("ID")
            self.FORMAT[id] = line_dict
        # Alternative allele
        if line.startswith("ALT"):
            line = line.removeprefix("ALT=")
            self.ALT.append(self._handle_xml_fmt(line))
        # assembly field
        if line.startswith("ass"):
            self.ASSEMBLY = line.removeprefix("assembly=")
        # contig field
        if line.startswith("con"):
            line = line.removeprefix("contig=")
            line_dict = self._handle_xml_fmt(line)
            id = line_dict.pop("ID")
            self.CONTIG[id] = line_dict
        # sample field
        if line.startswith("SAM"):
            line = line.removeprefix("SAMPLE=")
            self.SAMPLE.append(self._handle_xml_fmt(line))
        # pedigree
        if line.startswith("PEDIGREE="):
            line = line.removeprefix("PEDIGREE=")
            self.PEDIGREE.append(self._handle_xml_fmt(line))
        if line.startswith("pedigreeDB"):
            line = line.removeprefix("pedigreeDB=")
            self.PEDIGREE_DB = line
        pass

    def _handle_url_fmt(self, line) -> Tuple[str, str]:
        key, url = line.split("=")
        return key, url

    def _handle_xml_fmt(self, line) -> dict[str, str]:
        """
        Converts the XML-like format of VCF meta lines into a python Dictionary of strings
        example:  
        ```
        <
            ID=ID,
            Number=number,
            Type=type,
            Description="description",
            Source="source",
            Version="version"
        >
        ```
        is turned into  
        ```
        {  
            "ID": "ID",  
            "Number":"number",  
            "Type":"type",  
            "Description": "description",  
            "Source", "source",  
            "Version":"version  
        }
        ```  
        """
        print(line)
        fmt_line = line.strip("<>").split(",") #TODO: This causes commas within descriptions to be split, as well. Gotta change this...
        fmt_line = [x.split("=", maxsplit=1) for x in fmt_line]
        print(fmt_line)
        fmt_dict = {k:v.strip("\"") for k,v in fmt_line}
        return fmt_dict

    def __iter__(self):
        """loop over VCF records"""
        file = gzip.open(self.path, 'rt') if self.compressed else open(self.path, 'r')
        for line in file:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                yield Record(line)

class Record():
    def __init__(self, line:list):
        self.CHROM:str = line[0]
        self.POS:str = line[1]
        self.ID:str = line[2]
        self.REF:str = line[3]
        self.ALT:str = line[4]
        self.QUAL:str = line[5]
        self.FILTER:str = line[6]
        self.INFO:list = self._parse_INFO(line[7])
        try :
            self.FORMAT = self._parse_FORMAT(line[8])
            self.GT = self._parse_GT(line[9:])
        except IndexError:
            self.FORMAT = None
            self.GT = None

    def _parse_INFO(self, info_col:str) -> dict[str,list[str]]:
        """Parses INFO column and returns dictionary"""
        info_lst:list = info_col.split(";")
        info_kv_pairs:list = [ x.split("=") for x in info_lst ]
        return {k:v.split(",") for k,v in info_kv_pairs}
    
    def _parse_FORMAT(self, fmt_line) -> dict:
        """parses format column"""
        line = fmt_line.split(":")
        if type(line) == list:
            return {k:v for k,v in enumerate(line)}
        else:
            return {line:0}                             #handle cases where there is only one field
    
    #TODO redo this because genotype can be filled with multiple sub fields
    def _parse_GT(self, gt_cols:list[str]) -> list[Tuple[int, int, bool]]:
        """
        parses GT columns and returns List of Tuples with format (allele:int, allele:int, ... ,  phased?:bool)  
        missing genotypes are indicated as -1
        """
        # gt_array = []
        # phased = False
        # for col in gt_cols:
        #     col :str 
        #     phased = True if col[1] == "|" else False                           #check if phased
        #     alleles:list = col.split("|") if phased else col.split("/")         #split each genotype into alleles
        #     alleles = [ int(x) if x.isdigit() else -1 for x in alleles ]        #convert alleles into ints
        #     gt_array.append(alleles + [phased])
        # return gt_array
    
        return

if __name__ =="__main__":
    import time
    start_time = time.time()

    vcf = VCF("top50_000.vcf")
    # vcf = VCF("hprc-v2.0-mc-chm13.wave.year1.snps.vcf.gz", compressed=True)
    print(vcf.len)

    # for rec in vcf:
    #     print(rec.CHROM, rec.POS, rec.ID)

    end_time = time.time()
    print(f"took {end_time - start_time:.2f}s")