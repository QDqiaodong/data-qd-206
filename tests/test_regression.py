import unittest
from star_crossmatch import process
class Regression(unittest.TestCase):
 def test_wrap_and_invalid_threshold(self):
  d={'reference':[{'id':'r','ra':'359.9','dec':'0','magnitude':'1'}],'observations':[{'id':'o','ra':'0.1','dec':'0','sigma':'0.3'}],'threshold':'1'}
  r=process(d);self.assertEqual(r['matches'][0]['selected'],'r');self.assertAlmostEqual(float(r['matches'][0]['candidates'][0]['distance']),0.2,places=10)
  d['threshold']='0'
  with self.assertRaises(ValueError):process(d)
