import unittest

from hmnfusion.region import Region


class TestRegion(unittest.TestCase):
    """Test Region object"""

    # TODO: Add tests set_sequence_reference(), set_sequence_sample()
    def setUp(self):
        self.r_a = Region(
            chrom="chr7",
            position=100,
            start=75,
            end=125,
            interval=50,
            sequence_reference="ATCG",
            sequence_sample="ACCG",
            orientation="left",
        )
        self.d_a = dict(
            chrom="chr7",
            position=100,
            start=75,
            end=125,
            interval=50,
            sequence_reference="ATCG",
            sequence_sample="ACCG",
            orientation="left",
        )

    def test_init(self):
        """Test init()"""
        region = Region()
        self.assertEqual(region.chrom, "")
        self.assertEqual(region.position, 0)
        self.assertEqual(region.start, 0)
        self.assertEqual(region.end, 0)
        self.assertEqual(region.interval, 0)
        self.assertEqual(region.sequence_reference, "")
        self.assertEqual(region.sequence_sample, "")
        self.assertEqual(region.orientation, "undefined")

    def test_getters(self):
        """Test getters attributes"""
        self.assertEqual(self.r_a.chrom, "chr7")
        self.assertEqual(self.r_a.position, 100)
        self.assertEqual(self.r_a.start, 75)
        self.assertEqual(self.r_a.end, 125)
        self.assertEqual(self.r_a.interval, 50)
        self.assertEqual(self.r_a.sequence_reference, "ATCG")
        self.assertEqual(self.r_a.sequence_sample, "ACCG")
        self.assertEqual(self.r_a.orientation, "left")

    def test_setters(self):
        """Test setters attributes"""
        region = Region()
        region.chrom = "chr8"
        region.position = 10
        region.start = 11
        region.end = 12
        region.interval = 13
        region.sequence_reference = "A"
        region.sequence_sample = "T"
        region.orientation = "right"

        self.assertEqual(region.chrom, "chr8")
        self.assertEqual(region.position, 10)
        self.assertEqual(region.start, 11)
        self.assertEqual(region.end, 12)
        self.assertEqual(region.interval, 13)
        self.assertEqual(region.sequence_reference, "A")
        self.assertEqual(region.sequence_sample, "T")
        self.assertEqual(region.orientation, "right")

    def test_get_start(self):
        """Test get_start()"""
        self.assertEqual(self.r_a.get_start(), 75)
        self.r_a.start = 0
        self.assertEqual(self.r_a.start, 0)
        self.assertEqual(self.r_a.get_start(), 75)
        self.assertEqual(self.r_a.start, 75)
        # Interval not initialized
        self.r_a.start = 0
        self.r_a.interval = 0
        self.assertEqual(self.r_a.get_start(), 0)
        self.r_a.interval = 50
        # Position not initialized
        self.r_a.position = 0
        self.assertEqual(self.r_a.get_start(), 0)

    def test_get_end(self):
        """Test get_end()"""
        self.assertEqual(self.r_a.get_end(), 125)
        self.r_a.end = 0
        self.assertEqual(self.r_a.end, 0)
        self.assertEqual(self.r_a.get_end(), 125)
        self.assertEqual(self.r_a.end, 125)
        # Interval not initialized
        self.r_a.end = 0
        self.r_a.interval = 0
        self.assertEqual(self.r_a.get_end(), 0)
        self.r_a.interval = 50
        # Position not initialized
        self.r_a.position = 0
        self.assertEqual(self.r_a.get_end(), 0)

    def test_get_length(self):
        """Test get_length()"""
        self.assertEqual(self.r_a.get_length(), 50)
        self.r_a.interval = 0
        self.assertEqual(self.r_a.get_length(), 50)
        # Length negative
        self.r_a.start = 30
        self.r_a.end = 20
        self.assertEqual(self.r_a.get_length(), 0)
        # Length could not be calculated
        self.r_a.start = 30
        self.r_a.start = 0
        self.r_a.end = 0
        self.assertEqual(self.r_a.get_length(), 0)

    def test_format(self):
        """Test format()"""
        self.assertEqual(self.r_a.format(fusion=False), "chr7:75-125")
        self.r_a.end = 0
        self.assertEqual(self.r_a.format(fusion=False), "chr7:75-125")
        self.assertEqual(self.r_a.format(), "left chr7:100")

    def test_is_init(self):
        """Test is_init()"""
        region = Region()
        self.assertFalse(region.is_init())

        region = Region(chrom="chr7", position=100, orientation="left")
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

        self.assertNotEqual(a, self.r_a)
        a.chrom = "chr7"
        a.position = 100
        a.start = 75
        a.end = 125
        a.interval = 50
        a.sequence_reference = "ATCG"
        a.sequence_sample = "ACCG"
        a.orientation = "left"

        self.assertEqual(a, self.r_a)

    def test_to_dict(self):
        """Test to_dict()"""
        self.assertEqual(self.r_a.to_dict(), self.d_a)

    def test_from_dict(self):
        """Test from_dict()"""
        b = Region.from_dict(self.d_a)
        self.assertEqual(self.r_a, b)

    def test_check_region(self):
        """Test check_region()"""
        self.assertFalse(Region.check_region(":10"))
        self.assertFalse(Region.check_region("10-10"))
        self.assertTrue(Region.check_region("10:10"))
        self.assertTrue(Region.check_region("chr10:10"))

    def test_from_str(self):
        """Test from_str()"""
        r = Region.from_str("chr7:100")
        self.assertEqual(self.r_a.chrom, r.chrom)
        self.assertEqual(self.r_a.position, r.position)
        r = Region.from_str("chr7:75-125")
        self.assertEqual(self.r_a.chrom, r.chrom)
        self.assertEqual(self.r_a.start, r.start)
        self.assertEqual(self.r_a.end, r.end)


if __name__ == "__main__":
    unittest.main()
