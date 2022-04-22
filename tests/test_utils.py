import os
import tempfile
import unittest

from hmnfusion import utils
from main_test import Main_test


class TestUtils(Main_test):
    """Test functions in utils file"""

    def test_validate_name_sample(self):
        """Test validate_name_sample()"""
        self.assertTrue(utils.validate_name_sample("TestA"))
        self.assertTrue(utils.validate_name_sample("Test-A"))
        self.assertFalse(utils.validate_name_sample("Test A"))
        self.assertFalse(utils.validate_name_sample("Test	A"))

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


if __name__ == "__main__":
    unittest.main()
