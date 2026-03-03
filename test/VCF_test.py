from biofunc.vcf import VCF

class TestRecordGT:
    def test_meta_info(self):
        vcf = VCF("data/test1.vcf")
        assert vcf.fileformat == "VCFv4.2"
        # assert vcf.INFO == ""
        assert vcf.FILTER == {"ID":{"Description":"description"}}
        assert vcf.FORMAT == {"ID":{"Number":"number", "Type":"type", "Description":"description"}}
        assert vcf.assembly == "url"
        assert vcf.contig == {'ctg1': {'URL': 'ftp://somewhere.org/assembly.fa'}}
        assert vcf.PEDIGREE == {"Name_0":"G0-ID", "Name_1":"G1-ID", "Name_N":"GN-ID"}
        assert vcf.pedigreeDB == "URL"

        # TODO format ALT to also (possibly) split
        assert vcf.ALT == {'type': {'Description': 'description'}}
        # TODO format SAMPLE to split lists by semicolon ;
        assert vcf.SAMPLE == {'S_ID': {'Genomes': 'G1_ID;G2_ID;GK_ID', 'Mixture': 'N1;N2;NK', 'Description': 'S1;S2;SK'}}
        

    def test_additional_meta(self):
        vcf = VCF("data/test1.vcf")
        assert vcf.bcftools_viewVersion == "1.19+htslib-1.19"
        assert vcf.random_thing_here == "lalala"