import os
import unittest

import pandas as pd
import numpy as np
import example_data
from mss import alignment

data_path = os.path.dirname(os.path.join(example_data.__path__[0]))
batchpath = data_path + '/example_data/peakdata/mssoutput/'


class test_alignment(unittest.TestCase):
    def test_stack(self):
        """Tests stack functionlity"""
        # data_path = os.path.dirname(os.path.join(example_data.__path__[0]))
        # batchpath = data_path + '/batch/'

        a, b = alignment.stack(batchpath)
        assert len(b) == 2, "output size is wrong"
        assert type(b) == tuple, "output type is wrong"
        assert type(b[1].iat[6, 2]) == np.float32, "stacking seems off"

    def test_realignment(self):
        """Tests realignment functionality"""
        # data_path = os.path.dirname(os.path.join(example_data.__path__[0]))
        # batchpath = data_path + '/batch/'

        batch_name = 'test1'
        rt_error = 0.05
        MZ_error = 0.015
        file_type = 'csv'
        b = alignment.realignment(batchpath, batch_name, file_type,
                                  rt_error, MZ_error)
        assert len(b) > 0, "testing realingment length is off"
        assert type(b) == pd.core.frame.DataFrame, "output type is wrong"
        assert type(b.iat[1, 2]) == np.float32, "realignment seems to be off"
