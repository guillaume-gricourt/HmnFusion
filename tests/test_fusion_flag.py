import unittest

from hmnfusion.fusion_flag import FusionFlag


class TestFusionFlag(unittest.TestCase):
    """Class to test FusionFlag"""

    def test_is_interest(self):
        """Test is_interest()"""
        self.assertTrue(FusionFlag.is_interest(1))
        self.assertTrue(FusionFlag.is_interest(3))
        self.assertFalse(FusionFlag.is_interest(2))
        self.assertFalse(FusionFlag.is_interest(6))

    def test_is_consensus_primary(self):
        """Test is_consensus_primary()"""
        self.assertTrue(FusionFlag.is_consensus_primary(2))
        self.assertTrue(FusionFlag.is_consensus_primary(6))
        self.assertFalse(FusionFlag.is_consensus_primary(4))
        self.assertFalse(FusionFlag.is_consensus_primary(5))

    def test_is_consensus_secondary(self):
        """Test is_consensus_secondary()"""
        self.assertTrue(FusionFlag.is_consensus_secondary(4))
        self.assertTrue(FusionFlag.is_consensus_secondary(6))
        self.assertFalse(FusionFlag.is_consensus_secondary(3))
        self.assertFalse(FusionFlag.is_consensus_secondary(10))

    def test_is_consensus(self):
        """Test is_consensus()"""
        self.assertTrue(FusionFlag.is_consensus(4))
        self.assertTrue(FusionFlag.is_consensus(6))
        self.assertFalse(FusionFlag.is_consensus(8))
        self.assertFalse(FusionFlag.is_consensus(17))

    def test_is_genefuse(self):
        """Test is_genefuse()"""
        self.assertTrue(FusionFlag.is_genefuse(8))
        self.assertTrue(FusionFlag.is_genefuse(10))
        self.assertFalse(FusionFlag.is_genefuse(4))
        self.assertFalse(FusionFlag.is_genefuse(16))

    def test_is_hmnfusion(self):
        """Test is_hmnfusion()"""
        self.assertTrue(FusionFlag.is_hmnfusion(16))
        self.assertTrue(FusionFlag.is_hmnfusion(20))
        self.assertFalse(FusionFlag.is_hmnfusion(8))
        self.assertFalse(FusionFlag.is_hmnfusion(14))

    def test_is_lumpy(self):
        """Test is_lumpy()"""
        self.assertTrue(FusionFlag.is_lumpy(32))
        self.assertTrue(FusionFlag.is_lumpy(40))
        self.assertFalse(FusionFlag.is_lumpy(16))
        self.assertFalse(FusionFlag.is_lumpy(24))


if __name__ == "__main__":
    unittest.main()
