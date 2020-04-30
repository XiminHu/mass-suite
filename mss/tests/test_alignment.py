import os
import unittest

import pandas as pd 
import numpy as np
from tqdm import tqdm

#__import__("mass-suite")
import mss

data_path = os.path.join(mss.__path__[0], 'data')

class test_alignment(unittest.TestCase):
    
    def test_stack(self):
        """Tests stack functionlity"""
        
        batchpath = 'https://github.com/XiminHu/mass-suite/tree/master/example_data/batch/test1'
        
        assert type(listing) == np.array, "reference is wrong datatype"
        assert num_files == 2, "is this the correct test file?"
        
    def test_realignment(self):
        """Tests realignment functionality"""
        
        batchpath = 'https://github.com/XiminHu/mass-suite/tree/master/example_data/batch/'
        batch_name = 'test1'
        rt_error = 0.05
        MZ_error = 0.015
        
        assert type(mz_error) == int, "mz error is wrong type"
        assert type(RT_error) == int, "rt error is wrong type"
        assert type(alignment_df) == pandas.core.frame.DataFrame, "result is wrong datatype"
        