import filecmp
import os
import tempfile
import unittest

from hmnfusion import utils
from tests.main_test import Main_test


class TestUtils(Main_test):
    """Test functions in utils file"""

    # def abort()
    # def read_json()
    # def update_list()
    # def cmdline()
    # def check_bam_index()
    def test_check_fasta_index(self):
        """Test check_fasta_index()"""
        # is a fasta file
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as fod:
            fod.write(b">chr1\nATCG\n")
        self.assertTrue(utils.check_fasta_index(fod.name))
        self.assertTrue(utils.check_fasta_index(fod.name))
        os.remove(fod.name)
        # is not a fasta file
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as fod:
            fod.write(b"test")
        self.assertFalse(utils.check_fasta_index(fod.name))
        os.remove(fod.name)

    def test_find_executable(self):
        with self.assertRaises(utils.ExecutableNotFound):
            utils.find_executable(name="abcdefg")
        self.assertFalse(utils.find_executable(name="abcdefg", toraise=False))
        self.assertTrue(utils.find_executable(name="hmnfusion"))

    def test_validate_name_sample(self):
        """Test validate_name_sample()"""
        self.assertTrue(utils.validate_name_sample("TestA"))
        self.assertTrue(utils.validate_name_sample("Test-A"))
        self.assertFalse(utils.validate_name_sample("Test A"))
        self.assertFalse(utils.validate_name_sample("Test	A"))

    def test_bam_to_fastq_deflate(self):
        """Test bam_to_fastq()"""
        a, b = utils.bam_to_fastq(path=self.sample_m.bam, compress=0)
        self.assertTrue(a.endswith(".fastq"))
        self.assertTrue(b.endswith(".fastq"))
        self.assertTrue(filecmp.cmp(a, self.sample_m.fastq[0]))
        self.assertTrue(filecmp.cmp(b, self.sample_m.fastq[1]))
        os.remove(a)
        os.remove(b)

    def test_bam_to_fastq_inflate(self):
        """Test bam_to_fastq()"""
        a, b = utils.bam_to_fastq(path=self.sample_m.bam)
        self.assertTrue(a.endswith(".fastq.gz"))
        self.assertTrue(b.endswith(".fastq.gz"))

        self.assertTrue(Main_test.compare_file_gz(a, self.sample_m.fastq[0] + ".gz"))
        self.assertTrue(Main_test.compare_file_gz(b, self.sample_m.fastq[1] + ".gz"))
        os.remove(a)
        os.remove(b)


if __name__ == "__main__":
    unittest.main()
