import unittest

from hmnfusion.fusion import Fusion
from hmnfusion.region import Region

class TestRegion(unittest.TestCase):

	def setUp(self):
		# First
		self.fusion_1 = Fusion()

		# Second
		self.fusion_2 = Fusion()
		self.fusion_2_region_first = Region('chr9', 50, 'left')
		self.fusion_2_region_second = Region('chr22', 100, 'right')

		self.fusion_2.first = self.fusion_2_region_first
		self.fusion_2.second = self.fusion_2_region_second
		self.fusion_2.evidence = 20
		self.fusion_2.evidence_details = {'soft':10, 'split':50}
		self.fusion_2.depth = 100
		self.fusion_2.ident = 'LUM1'
		self.fusion_2.buildFrom = 'LUM0'
		self.fusion_2.isConsensus = True
		self.fusion_2.software = 'lumpy'

		self.fusion_2_dict = dict(software=['lumpy'], 
				first=self.fusion_2_region_first.to_dict(), 
				second=self.fusion_2_region_second.to_dict(), 
				evidence=20,
				evidence_details={'soft':10, 'split':50},
				ident='LUM1', 
				buildFrom=['LUM0'],
				isConsensus=True,
				depth=100)

		# Three
		self.fusion_3 = Fusion()
		self.fusion_3_region_first = Region('chr2', 500, 'left')
		self.fusion_3_region_second = Region('chr5', 2000, 'right')

		self.fusion_3.first = self.fusion_3_region_first
		self.fusion_3.second = self.fusion_3_region_second
		self.fusion_3.evidence = 50
		self.fusion_3.evidence_details = {'foo':100, 'bar':500}
		self.fusion_3.depth = 500
		self.fusion_3.ident = 'GEN1'
		self.fusion_3.buildFrom = 'GEN2'
		self.fusion_3.isConsensus = False
		self.fusion_3.software = 'genefuse'

		self.fusion_3_dict = dict(software=['genefuse'], 
				first=self.fusion_3_region_first.to_dict(), 
				second=self.fusion_3_region_second.to_dict(), 
				evidence=50,
				evidence_details={'foo':100, 'bar':500},
				ident='GEN1', 
				buildFrom=['GEN2'],
				isConsensus=False,
				depth=500)

	def test_init(self):
		self.assertEqual(self.fusion_1.first, Region())
		self.assertEqual(self.fusion_1.second, Region())
		self.assertEqual(self.fusion_1.evidence, 0)
		self.assertEqual(self.fusion_1.evidence_details, {})
		self.assertEqual(self.fusion_1.depth, 0)
		self.assertEqual(self.fusion_1.ident, '')
		self.assertEqual(self.fusion_1.buildFrom, set())
		self.assertEqual(self.fusion_1.isConsensus, False)
		self.assertEqual(self.fusion_1.software, set({'undefined'}))

	def test_getters(self):
		# Simple
		self.assertEqual(self.fusion_2.first, Region('chr9', 50, 'left'))
		self.assertEqual(self.fusion_2.second, Region('chr22', 100, 'right'))
		self.assertNotEqual(self.fusion_2.first, Region('chr22', 100, 'right'))

		self.assertEqual(self.fusion_3.evidence, 50)
		self.assertEqual(self.fusion_2.evidence_details, {'soft':10, 'split':50})
		self.assertEqual(self.fusion_2.depth, 100)
		self.assertEqual(self.fusion_2.ident, 'LUM1')
		self.assertTrue(self.fusion_2.isConsensus)

		# Complex
		self.assertIsInstance(self.fusion_2.buildFrom, set)
		self.assertEqual(self.fusion_2.buildFrom, set({'LUM0'}))

		self.assertIsInstance(self.fusion_2.software, set)
		self.assertEqual(self.fusion_2.software, set({'lumpy'}))
		
	def test_setters(self):
		# Simple
		self.fusion_2.first = self.fusion_3_region_first
		self.fusion_2.second = self.fusion_3_region_second

		self.assertEqual(self.fusion_2.first, self.fusion_3_region_first)
		self.assertEqual(self.fusion_2.second, self.fusion_3_region_second)
		self.assertNotEqual(self.fusion_2.first, self.fusion_2.second)

		self.fusion_2.evidence = self.fusion_3.evidence
		self.assertEqual(self.fusion_2.evidence, self.fusion_3.evidence)
		self.fusion_2.evidence_details = self.fusion_3.evidence_details
		self.assertEqual(self.fusion_2.evidence_details, self.fusion_3.evidence_details)
		self.fusion_2.depth = self.fusion_3.depth
		self.assertEqual(self.fusion_2.depth, self.fusion_3.depth)
		self.fusion_2.ident = self.fusion_3.ident
		self.assertEqual(self.fusion_2.ident, self.fusion_3.ident)
		self.fusion_2.isConsensus = self.fusion_3.isConsensus
		self.assertEqual(self.fusion_2.isConsensus, self.fusion_3.isConsensus)

		# Complex
		self.fusion_2.buildFrom = 'LUM1'
		self.assertEqual(self.fusion_2.buildFrom, set({'LUM0', 'LUM1'}))
		self.fusion_2.buildFrom = set({'LUM2', 'LUM3'})
		self.assertEqual(self.fusion_2.buildFrom, set({'LUM2', 'LUM3'}))
		self.fusion_2.buildFrom = ['LUM2', 'LUM3']
		self.assertEqual(self.fusion_2.buildFrom, set({'LUM2', 'LUM3'}))

		self.fusion_2.software = 'genefuse'
		self.assertEqual(self.fusion_2.software, set({'lumpy', 'genefuse'}))

		self.fusion_2.software = ['genefuse', 'other']
		self.assertEqual(self.fusion_2.software, set({'lumpy', 'genefuse', 'other'}))

	def test_is_near(self):
		raf = Region('chr5', 2000, 'right')
		ras = Region('chr6', 2000, 'left')
		rbf = Region('chr5', 2500, 'left')
		rbs = Region('chr10', 20, 'right')
		rbs2 = Region('chr6', 2200, 'left')

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

	def test_remove_from_name_cons(self):
		for buildfrom in ['CONS2', 'GEN5', 'LUM5']: 
			self.fusion_2.buildFrom = buildfrom
		self.assertEqual(self.fusion_2.buildFrom, set({'LUM0', 'CONS2', 'GEN5', 'LUM5'}))
		self.fusion_2.remove_fom_name_cons()
		self.assertEqual(self.fusion_2.buildFrom, set({'LUM0', 'GEN5', 'LUM5'}))
		
	def test_is_same_chrom(self):
		self.assertFalse(self.fusion_2.is_same_chrom())
		self.fusion_2.second = Region('chr9', 10, 'right')
		self.assertTrue(self.fusion_2.is_same_chrom())

	def test_set_region(self):
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
		self.assertEqual(self.fusion_2.first, self.fusion_2_region_first)
		self.assertEqual(self.fusion_2.second, self.fusion_2_region_second)

		self.fusion_2.swap_region()

		self.assertEqual(self.fusion_2.first, self.fusion_2_region_second)
		self.assertEqual(self.fusion_2.second, self.fusion_2_region_first)

	def test_to_dict(self):
		self.assertEqual(self.fusion_2.to_dict(), self.fusion_2_dict)
		self.assertEqual(self.fusion_3.to_dict(), self.fusion_3_dict)

	def test_from_dict(self):
		fusion_2 = Fusion.from_dict(self.fusion_2_dict)
		self.assertEqual(self.fusion_2, self.fusion_2)

		fusion_3 = Fusion.from_dict(self.fusion_3_dict)
		self.assertEqual(self.fusion_3, self.fusion_3)

	def test_greater_lether(self):
		self.assertTrue(self.fusion_2 < self.fusion_3)
		self.assertFalse(self.fusion_2 > self.fusion_3)

		self.fusion_2.evidence = 2000
		self.assertTrue(self.fusion_2 > self.fusion_3)
		self.assertFalse(self.fusion_2 < self.fusion_3)

	def test_equality(self):
		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		self.assertEqual(self.fusion_2, fusion_2_copy)

		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		fusion_2_copy.software = 'foo'
		self.assertNotEqual(self.fusion_2, fusion_2_copy)

		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		fusion_2_copy.ident = 'foo'
		self.assertNotEqual(self.fusion_2, fusion_2_copy)

		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		fusion_2_copy.swap_region()
		self.assertNotEqual(self.fusion_2, fusion_2_copy)

		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		fusion_2_copy.evidence = 0
		self.assertNotEqual(self.fusion_2, fusion_2_copy)

		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		fusion_2_copy.software = 'foo'
		fusion_2_copy.ident = 'foo'
		fusion_2_copy.swap_region()
		fusion_2_copy.evidence = 0
		self.assertNotEqual(self.fusion_2, fusion_2_copy)

		fusion_2_copy = Fusion.from_dict(self.fusion_2_dict)
		fusion_2_copy.evidence_details = {}
		self.assertNotEqual(self.fusion_2, fusion_2_copy)

if __name__ == '__main__':
	unittest.main()
