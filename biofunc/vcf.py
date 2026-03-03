import gzip
from pathlib import Path
import csv
from typing import Tuple, List

class VCF():
    def __init__(self, path:str, compressed:bool=False):
        self._path:str = path
        self._compressed:bool = compressed
        self.samples: list
        self._len : int = 0

        input_file = gzip.open(self._path, 'rt') if self._compressed else open(self._path, 'r')
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

    #region -------------------------- PROPERTIES -------------------------- 
    @property
    def fileformat(self):
        "VCF format version number"
        return getattr(self, "_fileformat", None)

    # TODO:
    @property
    def INFO(self) -> dict | None:
        ""
        return getattr(self, "_INFO", None)
    
    @property
    def FILTER(self) -> dict[str,dict] | None:
        """
        Filters that have been applied to the data  
        Returned in format :
        
        {  
            "ID":  {  "Description":  "desc"  },  
            "ID2": {  "Description":  "desc"  },  
            ...
        }
        """
        return getattr(self, "_FILTER", None)
    
    @property
    def FORMAT(self) -> dict[str, dict] | None:
        """
        Genotype field specifications  
        Returned in format :  

        {
            "ID" : {
                    "Number":  "number",
                    "Type": "type",
                    "Description": "desc"
                    },
            "ID2" : {
                    "Number":  "number",
                    "Type": "type",
                    "Description": "desc"
                    },
        }
        """
        return getattr(self, "_FORMAT", None)
    
    # TODO:
    @property
    def ALT(self):
        ""
        return getattr(self, "_ALT", None)
    
    @property
    def assembly(self) -> str | None:
        """
        URL field that specifies location of a fasta file containing breakpoint assemblies
        referenced in the VCF records for structural variants via the BKPTID INFO key.
        """
        return getattr(self, "_assembly", None)
    
    @property
    def contig(self) -> dict[str, dict] | None :
        """
        Description of contigs referred to in the VCF file
        """
        return getattr(self, "_contig", None)
    
    @property
    def SAMPLE(self):
        """
        Sample Genotype Mappings
        """
        return getattr(self, "_SAMPLE", None)
    
    @property
    def PEDIGREE(self) -> dict | None:
        """
        Relationships Between genoomes
        """
        return getattr(self, "_PEDIGREE", None)
    
    @property
    def pedigreeDB(self) -> str | None:
        """
        Link to database of relationships between genomes
        """
        return getattr(self, "_pedigreeDB", None)

    @property
    def len(self):
        """number of variants within the VCF file"""
        return self._len
    
    #endregion


    def _handle_header(self, cols:list):
        #first 8 columns should always be the same - CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
        #columns 9 is FORMAT, so col 10 and beyond must be sample data
        try:
            self.samples = cols[9:]
        except IndexError:
            self.samples = None


    #region -------------------------- META-INFO PARSERS -------------------------- 


    @staticmethod
    def _handle_url_fmt(line) -> Tuple[str, str]:
        key, url = line.split("=")
        return key, url

    @staticmethod
    def _handle_xml_fmt(line:str) -> dict[str, str]:
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
        fmt_line = line.strip("<>")
        row = next( csv.reader([fmt_line]) ) # convert string to list, skipping commas inside quotes
        kv_tuples = [x.partition("=")[::2] for x in row]
        return {k:v.strip('"') for k,v in kv_tuples} # convert list to dictionary

    @staticmethod
    def _handle_xml(self, line:str, attr_name):
        attr = getattr(self, attr_name, None)
        info = VCF._handle_xml_fmt(line)
        if attr == None:
            setattr(self, attr_name, info)
        else :
            attr = info

    @staticmethod
    def _handle_id_fmt(self, line:str, attr_name:str):
        """
        Handles xml-like formats with a known ID field  
        Returns Tuple of ID:str and Dictionary
        """
        fmt_dict = VCF._handle_xml_fmt(line)
        id = fmt_dict.pop("ID")
        attr = getattr(self, attr_name, None)
        if attr == None:
            setattr(self, attr_name, {id:fmt_dict})
        else:
            attr[id] = fmt_dict

    @staticmethod
    # TODO edit because there can be multiple IDs
    def _handle_alt(self, line:str, attr_name:str):
        """ALT fields ID's can be colon separated"""
        fmt_dict = VCF._handle_xml_fmt(line)
        id = fmt_dict.pop("ID")
        attr = getattr(self, attr_name, None)
        if attr == None:
            setattr(self, attr_name, {id:fmt_dict})
        else:
            attr[id] = fmt_dict
    
    @staticmethod
    def _handle_line(self, line:str, attr_name):
        attr = getattr(self, attr_name, None)
        if attr == None:
            setattr(self, attr_name, line)
        else:
            attr = line 

    # this dict returns functions to how VCF4.2 reserved metainfo should be handled
    metainfo_dict = {
        "fileformat":_handle_line,
        "INFO": _handle_id_fmt,
        "FILTER":_handle_id_fmt,
        "FORMAT":_handle_id_fmt,
        "ALT":_handle_alt,
        "assembly":_handle_line,
        "contig":_handle_id_fmt,
        "SAMPLE":_handle_id_fmt,
        "PEDIGREE":_handle_xml,
        "pedigreeDB":_handle_line,
        "default":_handle_line
    }

    def _handle_meta(self, line:str):
        pre_sep, sep, post_sep = line.strip().partition("=")
        # run the specific function from metainfo_dict
        try :
            meta_func = self.metainfo_dict[pre_sep]
            meta_func(self=self, line=post_sep, attr_name="_"+pre_sep)
        except KeyError :
            meta_func = self.metainfo_dict["default"]
            meta_func(self=self, line=post_sep, attr_name=pre_sep)
        
    #endregion


    def __iter__(self):
        """loop over VCF records"""
        file = gzip.open(self._path, 'rt') if self._compressed else open(self._path, 'r')
        for line in file:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                yield Record(line)

    def split(self, chunks=2) -> tuple:
        """
        Returns list of generators starting at different positions within the VCF file
        For easier use with multiprocessing
        """
        generators = Tuple()
        return generators


class Record():
    
    def __init__(self, line:list):
        #fixed fields
        self.CHROM:str = line[0]
        self.POS:str = int(line[1])
        self.ID:str = line[2]
        self.REF:str = line[3]
        self.ALT:List[str] = line[4].split(",") if "," in line[4] else line[4]
        self.QUAL:str = int(line[5])
        self.FILTER:str = line[6]

        #lazy load properties
        self._info_line = line[7]
        self._format_line = line[8]
        self._gt_line = line[9:]


    #region -------------------------- PROPERTIES -------------------------- 
    @property
    def INFO(self) -> dict[str, list[str]]:
        """Only parse INFO on request"""
        return self._parse_INFO(self._info_line)

    @property
    def FORMAT(self):
        """Only parse FORMAT field on request"""
        try :
            return self._parse_FORMAT(self._format_line)
        except :
            return None

    @property
    def GT_fields(self) -> None | dict :
        """Only parse GT fields on request"""
        if self.FORMAT == None:
            return None
        else :
            return self._parse_GT_field(self._gt_line)
        
    #endregion


    def _parse_INFO(self, info_col:str) -> dict[str,list[str]]:
        """Parses INFO column and returns dictionary"""
        info_lst:list = info_col.split(";")
        info_kv_pairs:list = [ x.split("=") for x in info_lst ]
        return {k:v.split(",") for k,v in info_kv_pairs}
    
    def _parse_FORMAT(self, fmt_line) -> dict:
        """parses format column"""
        line = fmt_line.split(":")
        if type(line) == list:
            return {gt_type:idx for idx,gt_type in enumerate(line)}
        else:
            return {line:0}     #handle cases where there is only one field
    
    
    #region -------------------------- GENOTYPE FIELD PARSING -------------------------- 
    
    def _parse_GT_field(self, gt_cols:list[str]) -> dict:
        if self.FORMAT == None:
            return None
        else :
            gt_fields = {}
            # loop over format indices
            for gt_keyword, idx in self.FORMAT.items():
                
                if len(self.FORMAT) == 1:
                    #if only 1 format field, gt_cols is a 1d array so we can pass it as the arg
                    gt_fields[gt_keyword] = self._handle_GT_field(gt_keyword, gt_cols)
                else:
                    # gt_cols is a 2d array so we need to extract only the field we want
                    gt_keyword_field = [x[idx] for x in gt_cols] 
                    gt_fields[gt_keyword] = self._handle_GT_field(gt_keyword, gt_keyword_field)
                    
        return gt_fields
    


    @staticmethod
    def _handle_gt(gt_cols:list[str]) -> list[list[int|bool]]:
        """
        Return alleles per sample + phasing as format:  
        [a1(str), a2(str), ..., phase(bool)]
        """
        return [x.split(x[1])+[True] if x[1] == "|" else x.split(x[1])+[False] for x in gt_cols]

    # TODO
    @staticmethod
    def _handle_ft(gt_cols:list[str]):
        pass
    
    @staticmethod
    def _handle_list_type(type=str):
        """
        Return function that converts list into selected type
        """
        def _handle_list(gt_cols:list[str]):
            if gt_cols[0] == ".":
                return None
            return [type(x) for x in gt_cols]
        return _handle_list
    
    
    # pre-define list functions once to be reused throughout
    _int_col_handler = _handle_list_type(int)
    _float_col_handler = _handle_list_type(float)
    _str_col_handler = _handle_list_type(str)
    # dictionary of functions to handle keyword columns
    _gt_keyword_handler = {
            "GT": _handle_gt,
            "DP": _int_col_handler,
            "FT": _handle_ft,
            "GL": _float_col_handler,
            "GLE": _str_col_handler, #double check if its list of strings bc docs give GLE=0:... as example
            "PL": _int_col_handler,
            "GP": _float_col_handler,
            "GQ": _int_col_handler,
            "HQ": _int_col_handler,
            "PS": _int_col_handler,
            "PQ": _int_col_handler,
            "EC": _int_col_handler,
            "MQ": _int_col_handler,
    }

    def _handle_GT_field(self, gt_keyword:str, gt_cols:list[str]):
        return self._gt_keyword_handler[gt_keyword](gt_cols)
    
    #endregion