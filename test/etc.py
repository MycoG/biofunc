from biofunc.vcf import VCF

a = VCF("data/test1.vcf")
# print("INFO", a.INFO)
# print("FILTER", a.FILTER)
# print("FORMAT", a.FORMAT)
print("ALT", a.ALT)
# print("assembly", a.assembly)
# print("contig", a.contig)
# print("sample", a.SAMPLE)
# print("PEDIGREE", a.PEDIGREE)
# print("pedigreeDB", a.pedigreeDB)
# print("bcftools_viewVersion", a.bcftools_viewVersion)