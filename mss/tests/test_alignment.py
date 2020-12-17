import os
import unittest

import pandas as pd
import numpy as np
import mss
from mss import align

data_path = os.path.dirname(os.path.join(mss.__path__[0]))
batchpath = data_path + '/example_data/'


class test_alignment(unittest.TestCase):
    def test_stack(self):
        """Tests stack functionlity"""
        test = align.mss_process(batchpath, None, thres_noise=20000,
        	                     enable_score=False)
        assert len(test) > 0, "no df read in the folder"
#        assert len(b) > 0, "output size is wrong"
#        assert type(b) == pd.core.frame.DataFrame, "output type is wrong"
#        assert b[b.columns[0]].dtypes == np.float32, "stacking seems off"
#        assert len(b.columns) > 0, 'no columns'

#    def test_realignment(self):
#        """Tests realignment functionality"""
#        batch_name = 'test1'
#        rt_error = 0.05
#        MZ_error = 0.015
#        file_type = '.csv'
#        b = align.mss_align(batchpath, batch_name, file_type,
#                                  rt_error, MZ_error)
#        assert len(b) > 0, "testing realingment length is off"
#        assert type(b) == pd.core.frame.DataFrame, "output type is wrong"
#        assert b[b.columns[0]].dtypes == np.float32, "wrong value type"
#        assert b.columns[0] == 'Average m/z', 'wrong column'
#        assert b.columns[1] == 'Average RT (min)', 'wrong column'
#        assert b.columns[2] == 'Average sn', 'wrong column'
#        assert b.columns[3] == 'Average score', 'wrong column'
