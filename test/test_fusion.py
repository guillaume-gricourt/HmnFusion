import unittest

from hmnfusion._version import __app_name__
from hmnfusion.evidence import Evidence
from hmnfusion.fusion import Fusion
from hmnfusion.region import Region


class TestFusion(unittest.TestCase):
    """Class to test Fusion"""

    @classmethod
    def setUp(self):
        """Initialize 3 Fusions object"""
        # First
        self.fusion_1 = Fusion()

        # Second
        self.evidence_2_d = dict(
            raw=50,
            split=50,
            mate=0,
            clipped=10,
            depth=100
        )
        self.evidence_2 = Evidence.from_dict(self.evidence_2_d)

        self.fusion_2 = Fusion("genefuse")
        self.fusion_2_region_first = Region("chr9", 50, "left")
        self.fusion_2_region_second = Region("chr22", 100, "right")

        self.fusion_2.first = self.fusion_2_region_first
        self.fusion_2.second = self.fusion_2_region_second
        self.fusion_2.evidence = self.evidence_2
        self.fusion_2.number = 1
        self.fusion_2.is_consensus = False

        self.fusion_2_dict = dict(
            first=self.fusion_2_region_first.to_dict(),
            second=self.fusion_2_region_second.to_dict(),
            evidence=self.evidence_2.to_dict(),
            number=1,
            is_consensus=False,
            software="genefuse"
        )

        # Three3
        self.evidence_3_d = dict(
            raw=50,
            split=100,
            mate=0,
            clipped=500,
            depth=500
        )
        self.evidence_3 = Evidence.from_dict(self.evidence_3_d)

        self.fusion_3 = Fusion()
        self.fusion_3_region_first = Region("chr2", 500, "left")
        self.fusion_3_region_second = Region("chr5", 2000, "right")

        self.fusion_3.first = self.fusion_3_region_first
        self.fusion_3.second = self.fusion_3_region_second
        self.fusion_3.evidence = self.evidence_3
        self.fusion_3.number = 2
        self.fusion_3.is_consensus = True
        self.fusion_3.software = "lumpy"

        self.fusion_3_dict = dict(
            first=self.fusion_3_region_first.to_dict(),
            second=self.fusion_3_region_second.to_dict(),
            evidence=self.evidence_3.to_dict(),
            number=2,
            is_consensus=True,
            software="lumpy"
        )

    def test_init(self):
        """Check initialize attributes values"""
        # Empty args.
        self.assertEqual(self.fusion_1.first, Region())
        self.assertEqual(self.fusion_1.second, Region())
        self.assertEqual(self.fusion_1.evidence, Evidence())
        self.assertEqual(self.fusion_1.number, 0)
        self.assertEqual(self.fusion_1.software, __app_name__)
        # With args.
        self.assertEqual(self.fusion_2.software, "genefuse")

    def test_getters(self):
        """Test getters attributes"""
        self.assertEqual(
            self.fusion_2.first,
            Region("chr9", 50, "left")
        )
        self.assertEqual(
            self.fusion_2.second,
            Region("chr22", 100, "right")
        )
        self.assertNotEqual(
            self.fusion_2.first,
            Region("chr22", 100, "right")
        )

        self.assertEqual(self.fusion_2.evidence, self.evidence_2)
        self.assertEqual(self.fusion_2.number, 1)
        self.assertFalse(self.fusion_2.is_consensus)
        self.assertEqual(self.fusion_2.software, "genefuse")

    def test_setters(self):
        """Test setters attributes"""
        self.fusion_2.first = self.fusion_3_region_first
        self.fusion_2.second = self.fusion_3_region_second

        self.assertEqual(self.fusion_2.first, self.fusion_3_region_first)
        self.assertEqual(self.fusion_2.second, self.fusion_3_region_second)
        self.assertNotEqual(self.fusion_2.first, self.fusion_2.second)

        self.fusion_2.evidence = self.fusion_3.evidence
        self.assertEqual(self.fusion_2.evidence, self.fusion_3.evidence)

        self.fusion_2.number = self.fusion_3.number
        self.assertEqual(self.fusion_2.number, self.fusion_3.number)

        self.fusion_2.is_consensus = True
        self.assertTrue(self.fusion_2.is_consensus)

        self.fusion_2.software = self.fusion_3.software
        self.assertEqual(self.fusion_2.software, self.fusion_3.software)

    def test_update(self):
        """Test update(self)"""
        self.fusion_2.update(self.fusion_3)
        self.assertEqual(self.fusion_2.first, self.fusion_3.first)
        self.assertEqual(self.fusion_2.second, self.fusion_3.second)
        self.assertEqual(self.fusion_2.evidence, self.fusion_3.evidence)
        self.assertEqual(self.fusion_2.number, 1)

    def test_get_name(self):
        """Test get_name(self)"""
        self.assertEqual(
            self.fusion_1.get_name(),
            __app_name__[:3].upper()+"_0"
        )
        self.assertEqual(self.fusion_2.get_name(), "GEN_1")
        self.assertEqual(self.fusion_3.get_name(), "HMN_2")

    def test_get_software(self):
        """Test get_software(self)"""
        self.assertEqual(self.fusion_1.get_software(), __app_name__)
        self.assertEqual(self.fusion_2.get_software(), "genefuse")
        self.assertEqual(self.fusion_3.get_software(), __app_name__)

    def get_is_consensus(self):
        """Test get is_consensus(self)"""
        self.assertTrue(self.fusion_1.is_consensus)
        self.assertFalse(self.fusion_2.is_consensus)
        self.assertFalse(self.fusion_3.is_consensus)

    def set_is_consensus(self):
        """Test set is_consensus(self)"""
        self.assertTrue(self.fusion_1.is_consensus)
        self.fusion_1.is_consensus = True
        self.assertTrue(self.fusion_1.is_consensus)
        self.assertEqual(self.fusion_1.software, __app_name__)

        self.assertFalse(self.fusion_2.is_consensus)
        self.fusion_2.is_consensus = True
        self.assertTrue(self.fusion_2.is_consensus)
        self.assertEqual(self.fusion_2.software, __app_name__)

    def test_is_near(self):
        """Test is_near()"""
        raf = Region("chr5", 2000, "right")
        ras = Region("chr6", 2000, "left")
        rbf = Region("chr5", 2500, "left")
        rbs = Region("chr10", 20, "right")
        rbs2 = Region("chr6", 2200, "left")

        # Check Chrom differents
        self.assertFalse(self.fusion_2.is_near(self.fusion_3, 10))
        # Check chrom first are same
        self.fusion_2.first = raf
        self.fusion_2.second = ras
        self.fusion_3.first = rbf
        self.fusion_3.second = rbs
        self.assertFalse(self.fusion_2.is_near(self.fusion_3, 10))
        # Check chrom equal but two breakpoints are too far
        self.fusion_3.second = rbs2
        self.assertFalse(self.fusion_2.is_near(self.fusion_3, 100))
        # Check chrom equal but one breakpoints is too far
        self.assertTrue(self.fusion_2.is_near(self.fusion_3, 250))
        # Check chrom equal but two breakpoints are in the target
        self.assertTrue(self.fusion_2.is_near(self.fusion_3, 501))

    def test_is_same_chrom(self):
        """Test is_same_chrom()"""
        self.assertFalse(self.fusion_2.is_same_chrom())
        self.fusion_2.second = Region("chr9", 10, "right")
        self.assertTrue(self.fusion_2.is_same_chrom())

    def test_set_region(self):
        """Test set_region()"""
        self.assertEqual(self.fusion_1.first, Region())
        self.assertEqual(self.fusion_1.second, Region())

        self.fusion_1.set_region(self.fusion_2_region_first)
        self.assertEqual(self.fusion_1.first, self.fusion_2_region_first)
        self.assertEqual(self.fusion_1.second, Region())

        self.fusion_1.set_region(self.fusion_2_region_first)
        self.assertEqual(self.fusion_1.first, self.fusion_1.second)

        self.fusion_1.set_region(self.fusion_2_region_second)
        self.assertEqual(self.fusion_1.first, self.fusion_2_region_first)
        self.assertEqual(self.fusion_1.second, self.fusion_2_region_second)

    def test_swap_region(self):
        """Test swap_region()"""
        self.assertEqual(self.fusion_2.first, self.fusion_2_region_first)
        self.assertEqual(self.fusion_2.second, self.fusion_2_region_second)

        self.fusion_2.swap_region()

        self.assertEqual(self.fusion_2.first, self.fusion_2_region_second)
        self.assertEqual(self.fusion_2.second, self.fusion_2_region_first)

    def test_to_dict(self):
        """Test to_dict()"""
        self.assertEqual(self.fusion_2.to_dict(), self.fusion_2_dict)
        self.assertEqual(self.fusion_3.to_dict(), self.fusion_3_dict)

    def test_from_dict(self):
        """Test from_dict()"""
        fusion_2 = Fusion.from_dict(self.fusion_2_dict)
        self.assertEqual(fusion_2, self.fusion_2)

        fusion_3 = Fusion.from_dict(self.fusion_3_dict)
        self.assertEqual(fusion_3, self.fusion_3)

    def test_greater_lether(self):
        """Test greather, lether metafunctions"""
        self.assertTrue(self.fusion_2 < self.fusion_3)
        self.assertFalse(self.fusion_2 > self.fusion_3)

        self.fusion_2.evidence = Evidence(2000)
        self.assertTrue(self.fusion_2 > self.fusion_3)
        self.assertFalse(self.fusion_2 < self.fusion_3)

    def test_equality(self):
        """Test equality metafunctions"""
        fusion_1b = Fusion()
        self.assertEqual(self.fusion_1, fusion_1b)
        self.assertNotEqual(self.fusion_1, self.fusion_2)


if __name__ == "__main__":
    unittest.main()
