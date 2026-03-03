from biofunc.vcf import VCF

a = VCF("test1.vcf")
print("INFO", a.INFO)
# print("FILTER", a.FILTER)
# print("FORMAT", a.FORMAT)
# print("ALT", a.ALT)
# print("assembly", a.assembly)