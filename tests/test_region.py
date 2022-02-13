import unittest

from hmnfusion.region import Region


class TestRegion(unittest.TestCase):
    """Test Region object"""

    def test_init(self):
        """Test init()"""
        region = Region()
        self.assertEqual(region.chrom, "")
        self.assertEqual(region.position, 0)
        self.assertEqual(region.orientation, "undefined")

    def test_getters(self):
        """Test getters attributes"""
        region = Region("chr7", 100, "left")
        self.assertEqual(region.chrom, "chr7")
        self.assertEqual(region.position, 100)
        self.assertEqual(region.orientation, "left")

    def test_setters(self):
        """Test setters attributes"""
        region = Region()
        region.chrom = "chr8"
        region.position = 10
        region.orientation = "right"

        self.assertEqual(region.chrom, "chr8")
        self.assertEqual(region.position, 10)
        self.assertEqual(region.orientation, "right")

    def test_is_init(self):
        """Test is_init()"""
        region = Region()
        self.assertFalse(region.is_init())

        region = Region("chr7", 100, "left")
        self.assertTrue(region.is_init())

        region.position = 0
        self.assertTrue(region.is_init())
        region.orientation = ""
        self.assertTrue(region.is_init())
        region.chrom = ""
        self.assertFalse(region.is_init())

    def test_orientation(self):
        """Test orientation attribute management"""
        region = Region()
        self.assertEqual(region.orientation, "undefined")
        region.orientation = "left"
        self.assertEqual(region.orientation, "left")
        region.orientation = "right"
        self.assertEqual(region.orientation, "right")
        region.orientation = "other"
        self.assertEqual(region.orientation, "undefined")

    def test_equal(self):
        """Test equal metafunction"""
        a = Region()
        b = Region("chr7", 100, "left")

        self.assertNotEqual(a, b)
        a.chrom = "chr7"
        a.position = 100
        a.orientation = "left"

        self.assertEqual(a, b)

    def test_to_dict(self):
        """Test to_dict()"""
        region = Region("chr7", 100, "left")
        data = dict(chrom="chr7", orientation="left", position=100)
        self.assertEqual(region.to_dict(), data)

    def test_from_dict(self):
        """Test from_dict()"""
        a = Region("chr7", 100, "left")
        data = dict(chrom="chr7", orientation="left", position=100)
        b = Region.from_dict(data)
        self.assertEqual(a, b)

    def test_check_region(self):
        """Test check_region()"""
        self.assertFalse(Region.check_region(":10"))
        self.assertFalse(Region.check_region("10-10"))
        self.assertTrue(Region.check_region("10:10"))
        self.assertTrue(Region.check_region("chr10:10"))


if __name__ == "__main__":
    unittest.main()
