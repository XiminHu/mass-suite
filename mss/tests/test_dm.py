import unittest
import os
import pandas as pd
import numpy as np
from mss import dm
import mss
test_path = os.path.join(mss.__path__[0], 'test')
data_path = os.path.join(test_path, 'data')
file = 'data/sample1114.csv'
file_path = os.path.join(data_path, file)

# Test update: how to write this type of test in chain
# output from one function is input for next function


class test_dm(unittest.TestCase):
    def test_data_prep(self):
        '''test data prep function'''
        d_ms = pd.read_csv(file_path)
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_test = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                              mz_range=[200, 800], area_thres=500,
                              simp_summary=False)
        assert len(d_test) == 5299, 'dataframe length error!'
        # More assert here
        return

    def test_ms_cluster(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        d_test = dm.ms_cluster(d_sample, ['SR520-Cal'], 'linear',
                               d_reduce=False, visual=False,
                               cluster_method='optics',
                               eps=0.6, min_samples=10)
        assert sum(d_test['label'] == np.nan) == 0, 'invalid label presents'
        assert len(list(set(d_test['label']))) == 16, 'incorrect clustering'
        assert d_test.shape[1]-d_sample.shape[1] == 1, 'column insert error'
        return

    def test_source_label(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        sourcelist = ['Coulter', 'Crescent', 'Miller',
                      'Swan', 'SR520-Cal-in-DI_1000mL']
        # Needs adjustment
        d_test = dm.source_label(d_sample, sourcelist,
                                 area_thres=50000, concat=True)
        assert d_test.shape[1] - d_sample.shape[1] == 11, 'column number error'
        assert d_test['source'].dtype == object, 'source column type error'
        return

    def test_source_report(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        sourcelist = ['Coulter', 'Crescent', 'Miller', 'Swan',
                      'SR520-Cal-in-DI_1000mL']
        # Needs adjustment
        d_label = dm.source_label(d_sample, sourcelist, area_thres=50000,
                                  concat=True)
        d_test = dm.source_report(d_label, ['Coulter', 'Crescent',
                                  'Miller', 'Swan', 'SR520-Cal-in-DI_1000mL'],
                                  ['Mix'], method='multiple',
                                  pa_thres=10000, CV_thres=2)
        assert d_test.shape == (27, 16), 'dataframe shape error'
        assert d_test['sample'].dtype == object, 'sample column type error'
        return
