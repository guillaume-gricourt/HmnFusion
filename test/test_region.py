import unittest

from hmnfusion.region import Region

class TestRegion(unittest.TestCase):

	def test_init(self):
		region = Region()
		self.assertEqual(region.chrom, '')
		self.assertEqual(region.position, 0)
		self.assertEqual(region.orientation, 'undefined')

	def test_getters(self):
		region = Region('chr7', 100, 'left')
		self.assertEqual(region.chrom, 'chr7')
		self.assertEqual(region.position, 100)
		self.assertEqual(region.orientation, 'left')

	def test_setters(self):
		region = Region()
		region.chrom = 'chr8'
		region.position = 10
		region.orientation = 'right'

		self.assertEqual(region.chrom, 'chr8')
		self.assertEqual(region.position, 10)
		self.assertEqual(region.orientation, 'right')

	def test_is_init(self):
		region = Region()
		self.assertFalse(region.is_init())

		region = Region('chr7', 100, 'left')
		self.assertTrue(region.is_init())

		region.position = 0
		self.assertTrue(region.is_init())
		region.orientation = ''
		self.assertTrue(region.is_init())
		region.chrom = ''
		self.assertFalse(region.is_init())

	def test_orientation(self):
		region = Region()
		self.assertEqual(region.orientation, 'undefined')
		region.orientation = 'left'
		self.assertEqual(region.orientation, 'left')
		region.orientation = 'right'
		self.assertEqual(region.orientation, 'right')
		region.orientation = 'other'
		self.assertEqual(region.orientation, 'undefined')

	def test_equal(self):
		a = Region()
		b = Region('chr7', 100, 'left')

		self.assertNotEqual(a, b)
		a.chrom = 'chr7'
		a.position = 100
		a.orientation = 'left'

		self.assertEqual(a, b)

	def test_to_dict(self):
		region = Region('chr7', 100, 'left')
		data = dict(chrom='chr7', orientation='left', position=100)
		self.assertEqual(region.to_dict(), data)

	def test_from_dict(self):
		a = Region('chr7', 100, 'left')
		data = dict(chrom='chr7', orientation='left', position=100)
		b = Region.from_dict(data)
		self.assertEqual(a,b)

if __name__ == '__main__':
	unittest.main()
