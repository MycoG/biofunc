from vcf_parser import VCF
from pysam import VariantFile
import cyvcf2
import time
from cProfile import Profile


if __name__ =="__main__":
    
    print("test vcf_parse")
    start_time = time.time()
    vcf = VCF("top50_000.vcf")
    with open("test_vcf_parse.bed", 'w') as f:
        for rec in vcf:
            f.write("\t".join([rec.CHROM, rec.POS, rec.ID]) + "\n")
    end_time = time.time()
    duration1 = end_time - start_time
    print(f"took {duration1:.2f}s")
    print(f"time for 20_000_000: {duration1*20_000_000/50_000:.3f}s")

    print("test pysam Variant")
    start_time = time.time()
    vcf = VariantFile("top50_000.vcf")
    with open("test_vcf_parse.pysam.bed", 'w') as f:
        for rec in vcf.fetch():
            f.write("\t".join([rec.chrom, str(rec.pos), rec.id]) + "\n")
    end_time = time.time()
    duration2 = end_time - start_time
    print(f"took {duration2:.2f}s")
    print(f"time for 20_000_000: {duration2*20_000_000/50_000:.3f}s")

    print("test cyVCF Variant")
    start_time = time.time()
    vcf = cyvcf2.VCF("top50_000.vcf")
    with open("test_vcf_parse.cyvcf2.bed", 'w') as f:
        for rec in vcf:
            f.write("\t".join([rec.CHROM, str(rec.POS), rec.ID]) + "\n")
    end_time = time.time()
    duration2 = end_time - start_time
    print(f"took {duration2:.2f}s")
    print(f"time for 20_000_000: {duration2*20_000_000/50_000:.3f}s")

    print(f"duration 1 2 difference : {duration1 - duration2:.2f}s")
    print(f"duration 1 2 difference 20_000_000 : {(duration1 - duration2)*20_000_000/50_000:.2f}s")


