import gzip
from pathlib import Path

class _Record():
    def __init__(self,
                 CHROM,
                 POS,
                 ID,
                 REF,
                 ALT,
                 QUAL,
                 FILTER,
                 INFO,
                 FORMAT,
                 GENOTYPES):
        self._chrom = CHROM
        self._pos = POS
        self._id = ID 
        self._ref = REF 
        self._alt = ALT 
        self._qual = QUAL
        self._filter = FILTER
        self._info = INFO
        self._format = FORMAT
        self._genotypes = GENOTYPES
    
    @property
    def chrom(self):
        return self._chrom
    
    @property
    def pos(self):
        return self._pos
    
    @property
    def id(self):
        return self._id
    
    @property
    def ref(self):
        return self._ref
    
    @property
    def alt(self):
        return self._alt
    
    @property
    def qual(self):
        return self._qual
    
    @property
    def filter(self):
        return self._filter
    
    @property
    def info(self):
        return self._info
    
    @property
    def format(self):
        return self._format
    
    @property
    def genotypes(self):
        return self._genotypes

class VCF():
    def __init__(self, path, compressed=False):
        self.path = Path(path) 
        self.compressed = compressed
        self.fileformat : str
        self.sampleids : list
        self.contigs = {}
        self.info = {}
        self.filter = {}

        if compressed:
            file = gzip.open(path, 'rt')
        else:
            file = open(path, 'r')

        #parse info lines at first
        for line in file:
            line = line.strip()
            if line.startswith("##"):
                self._parse_metainfo(line)
            elif line.startswith("#"):
                self._parse_header(line)
            else:
                break

        file.close()

    def _parse_field(self, field:str) -> dict:
        field = field.removeprefix("<").removesuffix(">")
        field_dict = {}
        items = field.split(",")
        for item in items:
            key, value = item.split("=")
            field_dict[key] = value
        return field_dict

    def _parse_metainfo(self, line:str):
        # remove the "##"
        line = line.removeprefix("##")
        key, value = line.split("=", 1)

        #match the key
        match key:
            case "INFO" :
                field_dict = self._parse_field(value)
                field_id = field_dict.pop("ID")
                self.info[field_id] = field_dict
            case "FILTER" :
                field_dict = self._parse_field(value)
                field_id = field_dict.pop("ID")
                self.filter[field_id] = field_dict
            case "FORMAT" :
                self.format = self._parse_field(value)
            case "ALT" :
                self.alt = self._parse_field(value)
            case "contig" :
                field_dict = self._parse_field(value)
                self.contigs[field_dict["ID"]] == field_dict["length"]
            case _ :
                pass
        return

    def _parse_header(self, line:str):
        line = line.removeprefix("#")
        try:
            self.sampleids = line.split("\t")[9:]
        except:
            print("No samples present...")
            self.sampleids = None

    def _parse_data_line(self, line:str) -> _Record:
        line = line.split("\t")
        return _Record(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9:])
    
    def _get_sample_idx(self, sample:str):
        return self.sampleids.index[sample]
    
    def iterrows(self):
        "generator to iterate over each record and return a _Record obj"
        if self.compressed:
            file = gzip.open(self.path, 'rt')
        else:
            file = open(self.path, 'r')

        for line in file:
            line = line.strip()
            if not line.startswith("#"):
                yield self._parse_data_line(line)
        return

    def to_bed(path, compressed=False):
        if compressed:
            gzip.open(path,'wt').write("").close()
            file = gzip.open(path, 'a')
        else:
            open(path, 'w').write("").close()
            file = open(path, 'r')

        #TODO: code for writing bed here

        file.close()
        return
