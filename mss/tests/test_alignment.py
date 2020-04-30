import os
import unittest

import pandas as pd 
import numpy as np
from tqdm import tqdm

import mss
from mss import alignment

data_path = os.path.join(mss.__path__[0], 'data')

class test_alignment(unittest.TestCase):
    
    def test_stack(self):
        """Tests stack functionlity"""
        
        batchpath = '/mnt/c/Users/nozom/desktop/mass-suite/example_data/batch/'
        a = alignment.stack(batchpath)
        
        assert len(a) == 2, "output size is wrong"
        assert type(a) == tuple, "output type is wrong"
        assert type(a[1].iat[6,2]) == np.float32, "stacking seems off"
        
    def test_realignment(self):
        """Tests realignment functionality"""
        
        batchpath = '/mnt/c/Users/nozom/desktop/mass-suite/example_data/batch/'
        batch_name = 'test1'
        rt_error = 0.05
        MZ_error = 0.015
        file_type = 'csv'
        
        # b = alignment.realignment(batchpath, batch_name, file_type, rt_error, MZ_error)

        # assert type(mz_error) == int, "mz error is wrong type"
        # assert type(RT_error) == int, "rt error is wrong type"
        # assert type(alignment_df) == pandas.core.frame.DataFrame, "result is wrong datatype"
        