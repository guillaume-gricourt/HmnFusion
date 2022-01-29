import os
import unittest


class Main_test(unittest.TestCase):
    dataset_path = os.path.join(os.path.dirname(__file__), "dataset")
    # HmnFusion.
    hmnfusion_path = os.path.join(dataset_path, "hmnfusion")
    test_a_ef = os.path.join(hmnfusion_path, "TEST-A.extractfusion.json")
    # Genefuse.
    genefuse_path = os.path.join(dataset_path, "genefuse")
    test_a_genefuse_html = os.path.join(genefuse_path, "TEST-A.html")
    # Lumpy.
    lumpy_path = os.path.join(dataset_path, "lumpy")
    test_a_lumpy = os.path.join(lumpy_path, "TEST-A.vcf")
