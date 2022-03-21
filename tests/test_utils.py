import unittest

from hmnfusion import utils
from main_test import Main_test


class TestUtils(Main_test):
    """Test functions in utils file"""

    def test_validate_name_sample(self):
        """Test validate_name_sample()"""
        with self.assertRaises(ValueError):
            utils.validate_name_sample("Test A")
        with self.assertRaises(ValueError):
            utils.validate_name_sample("Test	A")


if __name__ == "__main__":
    unittest.main()
