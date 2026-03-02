from biofunc.vcf import Record

record_line_gt = "chr10	2	>13041>13046_1	T	A,TCCCTCAC	60	.	AC=1,1;AF=0.0217391,0.0217391;AN=7;NS=232;LV=0;ORIGIN=chr10:1;LEN=1,7;TYPE=snp,ins	GT	.|0	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	2|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|1	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|."
record_line_gt = record_line_gt.strip().split("\t")

class TestRecordGT:
    
    def test_fixed_fields(self):
        rec = Record(record_line_gt)
        assert rec.CHROM == "chr10"
        assert rec.POS == 2
        assert rec.ID == ">13041>13046_1"
        assert rec.REF == "T"
        assert rec.ALT == ["A","TCCCTCAC"]
        assert rec.QUAL == 60
        assert rec.FILTER == "."
        
    def test_info(self):
        rec = Record(record_line_gt)
        info = rec.INFO
        assert info["AC"] == ["1","1"]
        assert info["AF"] == ["0.0217391","0.0217391"]
        assert info["AN"] == ["7"]
        assert info["NS"] == ["232"]
        assert info["LV"] == ["0"]
        assert info["ORIGIN"] == ["chr10:1"]
        assert info["LEN"] == ["1","7"]
        assert info["TYPE"] == ["snp","ins"]

        # TODO: convert str to ints and floats if required,
        # And if the item is single, do not put it in a list
        # assert info["AC"] == [1,1]
        # assert info["AF"] == [0.0217391,0.0217391]
        # assert info["AN"] == 7
        # assert info["NS"] == 232
        # assert info["LV"] == 0
        # assert info["ORIGIN"] == "chr10:1"
        # assert info["LEN"] == [1,7]
        # assert info["TYPE"] == ["snp","ins"]
        

    def test_format(self):
        rec = Record(record_line_gt)
        assert rec.FORMAT == {"GT":0}

    def test_gt(self):
        rec = Record(record_line_gt)
        gt = rec.GT_fields["GT"]
        assert gt[0] == [".","0", True]
        assert gt[-1] == [".", ".", True]

        # TODO: convert str to ints and "." to -1
        # assert gt[0] == [-1,0, True]
        # assert gt[-1] == [-1, -1, True]

class RecordTestNoGT:
    # test fixed fields
    # test format
    # test all gt_keywords
    pass