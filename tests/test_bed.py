import unittest

import pandas as pd
from hmnfusion.bed import Bed
from hmnfusion.region import Region
from main_test import Main_test


class TestBed(Main_test):
    """Test Bed object"""

    @classmethod
    def setUp(self):
        self.bed_a_df = pd.DataFrame(
            {"chrom" : ["1", "1", "2"], "start" : [5, 15, 5], "end" : [10, 25, 10]}
        )

    def test_init(self):
        """Test init()"""
        bed = Bed()
        self.assertTrue(bed.df.empty)

    def test_getters(self):
        """Test getters attributes"""
        bed = Bed()
        self.assertTrue(bed.df.empty)

    def test_setters(self):
        """Test setters attributes"""
        bed = Bed()
        bed.df = self.bed_a_df

        self.assertTrue(bed.df.compare(self.bed_a_df, keep_equal=True).empty)

    def test_from_bed(self):
        """Test from_bed()"""
        bed = Bed.from_bed(self.bed_a_path)
        self.assertTrue(bed.df.compare(self.bed_a_df).empty)

    def test_select_bed(self):
        """Test select_bed()"""
        # Init.
        bed = Bed()
        bed.df = self.bed_a_df

        # Test 1
        region = Region("1", 8, "left")
        sel = bed.df.apply(
            Bed.select_bed, axis=1, args=(region,)
        )
        self.assertEqual(sel.sum(), 1)

        # Test 2
        region = Region("1", 2, "left")
        sel = bed.df.apply(
            Bed.select_bed, axis=1, args=(region,)
        )
        self.assertEqual(sel.sum(), 0)


if __name__ == "__main__":
    unittest.main()
