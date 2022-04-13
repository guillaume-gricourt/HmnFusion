import os
import subprocess
import unittest
from collections import namedtuple


class Main_test(unittest.TestCase):
    Sample = namedtuple(
        "Sample",
        [
            "name",
            "bam",
            "sam",
            "genefuse_json",
            "genefuse_html",
            "lumpy",
            "extractfusion",
            "quantification",
            "region",
        ],
    )

    dataset_path = os.path.join(os.path.dirname(__file__), "dataset")
    # Sample.
    ref_path = os.path.join(dataset_path, "ref")
    bam_path = os.path.join(dataset_path, "bam")
    vcf_path = os.path.join(dataset_path, "vcf")
    hmnfusion_path = os.path.join(dataset_path, "hmnfusion")
    genefuse_path = os.path.join(dataset_path, "genefuse")
    lumpy_path = os.path.join(dataset_path, "lumpy")

    sample_a = Sample(
        name="TEST-A",
        bam=os.path.join(bam_path, "TEST-A.bam"),
        sam=os.path.join(bam_path, "TEST-A.sam"),
        genefuse_json="",
        genefuse_html=os.path.join(genefuse_path, "TEST-A.html"),
        lumpy=os.path.join(lumpy_path, "TEST-A.vcf"),
        extractfusion=os.path.join(hmnfusion_path, "TEST-A.extractfusion.json"),
        quantification="",
        region="",
    )
    sample_m = Sample(
        name="TEST-M",
        bam=os.path.join(bam_path, "TEST-M.bam"),
        sam="",
        genefuse_json=os.path.join(genefuse_path, "TEST-M.json"),
        genefuse_html=os.path.join(genefuse_path, "TEST-M.html"),
        lumpy=os.path.join(lumpy_path, "TEST-M.vcf"),
        extractfusion=os.path.join(hmnfusion_path, "TEST-M.extractfusion.json"),
        quantification=os.path.join(hmnfusion_path, "TEST-M.quantification.vcf"),
        region="",
    )
    sample_p = Sample(
        name="TEST-P",
        bam=os.path.join(bam_path, "TEST-P.bam"),
        sam=os.path.join(bam_path, "TEST-P.sam"),
        genefuse_json=os.path.join(genefuse_path, "TEST-P.json"),
        genefuse_html=os.path.join(genefuse_path, "TEST-P.html"),
        lumpy=os.path.join(lumpy_path, "TEST-P.vcf"),
        extractfusion=os.path.join(hmnfusion_path, "TEST-P.extractfusion.json"),
        quantification=os.path.join(hmnfusion_path, "TEST-P.quantification.vcf"),
        region="chr22:23632550",
    )

    # mmej-deletion
    ref_mmej = os.path.join(ref_path, "mmej.fasta")
    n1_vcf = os.path.join(vcf_path, "N1.vcf")
    n2_vcf = os.path.join(vcf_path, "N2.vcf")
    p1_vcf = os.path.join(vcf_path, "P1.vcf")
    p2_vcf = os.path.join(vcf_path, "P2.vcf")
    mmej_deletion_n1 = os.path.join(hmnfusion_path, "N1.mmej-deletion.xlsx")
    mmej_deletion_n1n2 = os.path.join(hmnfusion_path, "N1N2.mmej-deletion.xlsx")
    mmej_deletion_n1p1p2 = os.path.join(hmnfusion_path, "N1P1P2.mmej-deletion.xlsx")

    # Auxiliary.
    bed_bcr_path = os.path.join(dataset_path, "bed", "bcr.bed")
    bed_a_path = os.path.join(dataset_path, "bed", "a.bed")

    @classmethod
    def launch(cls, args):
        if isinstance(args, str):
            args = args.split()
        ret = subprocess.run(args, capture_output=True, encoding="utf8")
        return ret
