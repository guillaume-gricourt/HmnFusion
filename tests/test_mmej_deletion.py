import tempfile
import unittest

import pysam
import numpy as np
import pandas as pd
from hmnfusion import mmej_deletion
from main_test import Main_test


class TestMmejDeletionMain(Main_test):

    @classmethod
    def load_records(cls, path: str):
        vcf_in = pysam.VariantFile(path)
        return [x for x in vcf_in.fetch()]

    def setUp(self):
        # Value
        self.value_0 = mmej_deletion.Value()
        self.value_1 = mmej_deletion.Value(
            id="86ad494080bc9c322a639d3de922e958",
            contig="chr1",
            start=5,
            deletion="TGAGGC",
        )
        self.value_2 = mmej_deletion.Value(
            id="927f1d86b6d899d163efdb245b9aca67",
            contig="chr19",
            start=5,
            deletion="TGA"
        )
        self.value_1_df = pd.DataFrame(
            {
                "contig": "chr1",
                "start": 5,
                "deletion": "TGAGGC",
                "sequence": "TGAGGC",
                "conclusion": "alignment ambiguous"
            },
            index=["86ad494080bc9c322a639d3de922e958"],
        )
        self.value_2_df = pd.DataFrame(
            {
                "contig": "chr19",
                "start": 5,
                "deletion": "TGA",
                "sequence": "",
                "conclusion": "no clear signature"
            },
            index=["927f1d86b6d899d163efdb245b9aca67"],
        )

        self.values_unit_one = TestMmejDeletionMain.load_records(path=self.u1_vcf)

        # MmejDeletion
        self.mmej_deletion_u0 = mmej_deletion.MmejDeletion(name="sample0", values=[])
        self.mmej_deletion_u1 = mmej_deletion.MmejDeletion(
            name="sample1",
            values=[self.value_1],
        )
        self.mmej_deletion_u2_s1 = mmej_deletion.MmejDeletion(
            name="sample1",
            values=[self.value_1, self.value_2],
        )
        self.mmej_deletion_u2_s2 = mmej_deletion.MmejDeletion(
            name="sample2",
            values=[self.value_1],
        )
        self.mmej_deletion_u2_df = pd.concat([self.value_1_df, self.value_2_df])
        self.mmej_deletion_u2_df["sample1"] = ["o", "o"]
        self.mmej_deletion_u2_df["sample2"] = ["o", pd.NA]
        self.mmej_deletion_u2_df_xlsx = self.mmej_deletion_u2_df.replace({pd.NA: np.nan, "": np.nan})
        self.mmej_deletion_u2_df_xlsx.reset_index(inplace=True, drop=True)
        self.mmej_deletion_empty_df = pd.DataFrame(columns=["contig", "start", "deletion", "sequence", "conclusion", "N1"])
        self.mmej_deletion_empty_df_xlsx = pd.DataFrame({"Unnamed: 0": "no deletion found", "contig": np.nan, "start": np.nan, "deletion": np.nan, "sequence": np.nan, "conclusion": np.nan, "N1": np.nan}, index=[0])

class TestConclude(Main_test):
    """Test Conclude object"""

    def test_attribute(self):
        """Test attribute number"""
        attrs = [x for x in dir(mmej_deletion.Conclude) if not x.startswith("__")]
        self.assertEqual(len(attrs), 4)


class TestValue(TestMmejDeletionMain):
    """Test Value object"""

    def test_getters(self):
        """Test getters attributes"""
        self.assertEqual(self.value_1.id, "86ad494080bc9c322a639d3de922e958")
        self.assertEqual(self.value_1.contig, "chr1")
        self.assertEqual(self.value_1.start, 5)
        self.assertEqual(self.value_1.deletion, "TGAGGC")

    def test_setters(self):
        """Test setters attributes"""
        self.value_0.id = self.value_1.id
        self.value_0.contig = self.value_1.contig
        self.value_0.start = self.value_1.start
        self.value_0.deletion = self.value_1.deletion
        self.value_0.sequence = self.value_1.sequence
        self.assertEqual(self.value_0.id, "86ad494080bc9c322a639d3de922e958")
        self.assertEqual(self.value_0.contig, "chr1")
        self.assertEqual(self.value_0.start, 5)
        self.assertEqual(self.value_0.deletion, "TGAGGC")

    def test_get_conclusion(self):
        """Test get_conclusion()"""
        self.value_1.sequence = "ATCG"
        self.value_1.deletion = "ATCG"
        self.assertEqual(
            self.value_1.get_conclusion(),
            mmej_deletion.Conclude.AMBIGUOUS,
        )
        self.value_1.sequence = "ATC"
        self.assertEqual(
            self.value_1.get_conclusion(),
            mmej_deletion.Conclude.UNCLEAR,
        )
        self.value_1.sequence = "ATCGGC"
        self.assertEqual(
            self.value_1.get_conclusion(),
            mmej_deletion.Conclude.VALID,
        )
        self.value_1.deletion = "A"
        self.assertEqual(
            self.value_1.get_conclusion(),
            mmej_deletion.Conclude.UNINITIALIZED,
        )

    def test_set_sequence(self):
        """Test set_sequence()"""
        self.value_1.set_sequence(path=self.ref_mmej)
        self.assertEqual(self.value_1.sequence, "TGAGGC")

    def test_from_record(self):
        """Test from_record()"""
        rec = mmej_deletion.Value.from_record(self.values_unit_one[0])
        self.assertEqual(rec, self.value_1)

    def test_to_dataframe(self):
        """Test to_dataframe()"""
        self.value_1.set_sequence(path=self.ref_mmej)
        self.assertTrue(
            self.value_1.to_dataframe().equals(self.value_1_df)
        )

    def test_to_region(self):
        """Test to_region()"""
        self.assertEqual(self.value_1.to_region(), "chr1:5-17")


class TestMmejDeletion(TestMmejDeletionMain):
    """Test MmmejDeletion object"""

    def test_getters(self):
        """Test getters attributes"""
        self.assertEqual(self.mmej_deletion_u1.name, "sample1")
        self.assertEqual(self.mmej_deletion_u1.values, [self.value_1])

    def test_setters(self):
        """Test setters attributes"""
        self.assertEqual(self.mmej_deletion_u0.name, "sample0")
        self.assertEqual(self.mmej_deletion_u0.values, [])
        self.mmej_deletion_u0.name = self.mmej_deletion_u1.name
        self.mmej_deletion_u0.values = self.mmej_deletion_u1.values
        self.assertEqual(self.mmej_deletion_u1.name, "sample1")
        self.assertEqual(self.mmej_deletion_u1.values, [self.value_1])

    def test_empty(self):
        """Test empty property"""
        self.assertTrue(self.mmej_deletion_u0.empty)
        self.assertFalse(self.mmej_deletion_u1.empty)

    def test_build_empty_dataframe(self):
        """Test build_empty_dataframe"""
        self.assertTrue(mmej_deletion.MmejDeletion.build_empty_dataframe(name="test").equals(pd.DataFrame(columns=["contig", "start", "deletion", "sequence", "conclusion", "test"])))

    def test_get_value_ids(self):
        """Test get_value_ids()"""
        self.assertEqual(self.mmej_deletion_u0.get_value_ids(), [])
        self.assertEqual(
                self.mmej_deletion_u2_s1.get_value_ids(),
                ["86ad494080bc9c322a639d3de922e958", "927f1d86b6d899d163efdb245b9aca67"],
        )

    def test_set_value_sequence(self):
        """Test set_value_sequence()"""
        self.mmej_deletion_u0.set_value_sequence(path=self.ref_mmej)
        self.assertEqual(self.mmej_deletion_u0.values, [])

        self.assertEqual(self.mmej_deletion_u1.values[0].sequence, "")
        self.mmej_deletion_u1.set_value_sequence(path=self.ref_mmej)
        self.assertEqual(self.mmej_deletion_u1.values[0].sequence, "TGAGGC")

    def test_from_vcf(self):
        """Test from_vcf()"""
        dels = mmej_deletion.MmejDeletion.from_vcf(path=self.n1_vcf)
        self.assertEqual(dels, [mmej_deletion.MmejDeletion(name="N1", values=[])])

        dels = mmej_deletion.MmejDeletion.from_vcf(path=self.u2_vcf)
        self.assertEqual(
            dels,
            [self.mmej_deletion_u2_s1, self.mmej_deletion_u2_s2],
        )

    def test_to_dataframe(self):
        """Test to_dataframe()"""
        # Empty
        mmej_deletions = mmej_deletion.MmejDeletion.from_vcf(path=self.n1_vcf)
        for m in mmej_deletions:
            m.set_value_sequence(path=self.ref_mmej)
        df = mmej_deletion.MmejDeletion.to_dataframe(mmej_deletions=mmej_deletions)
        self.assertTrue(self.mmej_deletion_empty_df.equals(df))
        # Filled
        mmej_deletions = mmej_deletion.MmejDeletion.from_vcf(path=self.u2_vcf)
        for m in mmej_deletions:
            m.set_value_sequence(path=self.ref_mmej)
        df = mmej_deletion.MmejDeletion.to_dataframe(mmej_deletions=mmej_deletions)
        self.assertTrue(self.mmej_deletion_u2_df.equals(df))

    def test_to_excel(self):
        """Test to_excel()"""
        # Empty
        mmej_deletions = mmej_deletion.MmejDeletion.from_vcf(path=self.n1_vcf)
        for m in mmej_deletions:
            m.set_value_sequence(path=self.ref_mmej)
        with tempfile.NamedTemporaryFile(suffix=".xlsx") as fod:
            mmej_deletion.MmejDeletion.to_excel(path=fod.name, mmej_deletions=mmej_deletions)
            df = pd.read_excel(fod.name)
        self.assertTrue(self.mmej_deletion_empty_df_xlsx.equals(df))
        # Filled
        mmej_deletions = mmej_deletion.MmejDeletion.from_vcf(path=self.u2_vcf)
        for m in mmej_deletions:
            m.set_value_sequence(path=self.ref_mmej)
        with tempfile.NamedTemporaryFile(suffix=".xlsx") as fod:
            mmej_deletion.MmejDeletion.to_excel(path=fod.name, mmej_deletions=mmej_deletions)
            df = pd.read_excel(fod.name)
        self.assertTrue(self.mmej_deletion_u2_df_xlsx.equals(df))


if __name__ == "__main__":
    unittest.main()
