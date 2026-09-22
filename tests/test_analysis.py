import sys,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from analyze import validate_counts,load_data,enrich
class AnalysisTests(unittest.TestCase):
    def test_invalid_counts(self):
        for value in [-1,.5,np.nan,np.inf]:
            with self.assertRaises(ValueError): validate_counts(pd.DataFrame([[value]]))
    def test_duplicates(self):
        with self.assertRaises(ValueError): validate_counts(pd.DataFrame([[1,2]],columns=['A','A']))
    def test_paired_source_cohort(self):
        counts,meta=load_data()
        self.assertEqual(counts.shape,(26,26363))
        self.assertEqual(meta.donor.nunique(),13)
        self.assertTrue(counts.index.equals(meta.index))
        self.assertEqual(meta.gsm.nunique(),26)
    def test_enrichment_uses_measured_background(self):
        genes=[f'g{i}' for i in range(100)]
        de=pd.DataFrame({'padj':[.001]*10+[.9]*90,'log2FoldChange':[2]*10+[0]*90},index=genes)
        go=enrich(de,{'term':set(genes[:20])|{'unmeasured'}})
        row=go[go.direction=='Higher nasal'].iloc[0]
        self.assertEqual((row.M,row.N,row.K,row.k),(100,10,20,10))
        expected=fisher_exact([[10,0],[10,80]],alternative='greater').pvalue
        self.assertAlmostEqual(row.pvalue,expected,places=12)
        self.assertAlmostEqual(row.padj,min(1,2*expected),places=12)
    def test_empty_hit_set(self):
        genes=[f'g{i}' for i in range(20)]
        de=pd.DataFrame({'padj':[1.]*20,'log2FoldChange':[0.]*20},index=genes)
        go=enrich(de,{'term':set(genes)})
        self.assertTrue((go.padj==1).all())
if __name__=='__main__': unittest.main()
