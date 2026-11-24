import json,unittest
from pathlib import Path
from star_crossmatch import process
class T(unittest.TestCase):
 def setUp(self):self.d=json.loads(Path('fixtures/input.json').read_text())
 def test_match_and_unmatched(self):
  r=process(self.d);self.assertEqual(r['summary']['matched'],1);self.assertEqual(r['matches'][0]['selected'],'s1');self.assertEqual(r['matches'][1]['status'],'unmatched')
 def test_ambiguous(self):
  self.d['reference'][1]['ra']='10.2';self.d['observations'][0]['sigma']='1';r=process(self.d);self.assertEqual(r['matches'][0]['status'],'ambiguous')
 def test_duplicate_bad(self):
  self.d['observations'].append(self.d['observations'][0])
  with self.assertRaises(ValueError):process(self.d)
 def test_threshold(self):
  self.d['threshold']='0.1';self.assertEqual(process(self.d)['summary']['matched'],0)
 def test_missing_field(self):
  self.d['reference'][0].pop('magnitude')
  with self.assertRaises(ValueError):process(self.d)
 def test_tie_determinism(self):
  self.d['reference'][1]['ra']='9.9';self.d['observations'][0]['sigma']='1';self.assertEqual(process(self.d)['matches'][0]['candidates'][0]['reference'],'s1')
