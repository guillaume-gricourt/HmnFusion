import unittest

from hmnfusion.evidence import Evidence

class TestEvidence(unittest.TestCase):
    """Test Evidence object"""
    @classmethod
    def setUpClass(cls):
        """Setup objects used in test"""
        cls._e1a = Evidence()
        cls._e1a.raw = 10
        cls._e1a.split = 10        
        cls._e1a.mate = 10
        cls._e1a.clipped = 10
        cls._e1a.depth = 10

        cls._e1b = Evidence()
        cls._e1b.raw = 10
        cls._e1b.split = 10        
        cls._e1b.mate = 10
        cls._e1b.clipped = 10
        cls._e1b.depth = 10

        cls._e2 = Evidence(15)
        cls._e2.split = 40        
        cls._e2.mate = 0
        cls._e2.clipped = -100
        cls._e2.depth = -150
        cls._e2d = dict(raw=15, split=40, mate=0, clipped=100, depth=150)     

        cls._e3 = Evidence()
        cls._e4 = Evidence(20)
                    
	def test_init_empty(self):
		"""Test init() empty"""
		self.assertEqual(self._e3.raw,0)
		self.assertEqual(self._e3.split,0)
        self.assertEqual(self._e3.mate,0)
        self.assertEqual(self._e3.clipped,0)
        self.assertEqual(self._e3.depth,0)
        
	def test_init_filled(self):
		"""Test init() filled"""
		self.assertEqual(self._e4.raw,20)
		self.assertEqual(self._e4.split,0)
        self.assertEqual(self._e4.mate,0)
        self.assertEqual(self._e4.clipped,0)
        self.assertEqual(self._e4.depth,0)
        
	def test_getters(self):
		"""Test getters attributes"""        
		self.assertEqual(self._e1b.raw, 10)
		self.assertEqual(self._e1b.split, 10)
		self.assertEqual(self._e1b.mate, 10)
		self.assertEqual(self._e1b.clipped, 10)
		self.assertEqual(self._e1b.depth, 10)
		
		self.assertEqual(self._e2.raw, 15)
		self.assertEqual(self._e2.split, 40)
		self.assertEqual(self._e2.mate, 0)
		self.assertEqual(self._e2.clipped, 100)										
		self.assertEqual(self._e2.depth, 150)
		
	def test_setters(self):
		"""Test setters attributes"""
        cls._e3.raw = 0
        cls._e3.split = 1      
        cls._e3.mate = 2
        cls._e3.clipped = 3
        cls._e3.depth = 4

		self.assertEqual(self._e3.raw, 0)
		self.assertEqual(self._e3.split, 1)
		self.assertEqual(self._e3.mate, 2)
		self.assertEqual(self._e3.clipped, 3)
		self.assertEqual(self._e3.depth, 4)												

	def test_set_number(self):
		"""Test set_number"""
		self.assertEqual(Evidence.set_number("10"), 10)
		self.assertEqual(Evidence.set_number("-10"), 10)
		self.assertEqual(Evidence.set_number(15), 15)
		self.assertEqual(Evidence.set_number(-300), -300)
		self.assertEqual(Evidence.set_number("0.02"), 0)
		self.assertEqual(Evidence.set_number("1.5"), 1)
												
	def test_equal(self):
	    """Test equal"""
		self.assertEqual(self._e1a, self._e1b)
		self.assertNotEqual(self._e1a, self._e2)

    def test_get_vaf(self):
        """Test vaf fuction"""
        self.assertEqual(self._e1a.get_vaf(), '3.00')
        self.assertEqual(self._e1a.get_vaf(float), 3)
        self.assertEqual(self._e1a.get_vaf(int), 3)
        
        self.assertEqual(self._e2.get_vaf(), '0.93')
        self.assertEqual(self._e2.get_vaf(float), 0.93)
        self.assertEqual(self._e2.get_vaf(int), 0.93)
                        
	def test_to_dict(self):
		"""Test to_dict()"""     
		self.assertEqual(self._e2.to_dict(), self._e2d)

	def test_from_dict(self):
		"""Test from_dict()"""
		self.assertEqual(self._e2, Evidence.from_dict(self._e2d))

if __name__ == '__main__':
	unittest.main()
