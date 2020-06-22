import unittest
import os
import pandas as pd
import numpy as np
from mss import dm
import mss
import sklearn
import example_data
test_path = os.path.join(mss.__path__[0], 'tests')
data_path = os.path.join(test_path, 'data')
file = 'sample1114.csv'
file_path = os.path.join(data_path, file)
test_path2 = os.path.join(example_data.__path__[0], 'clustering')
file_path2 = os.path.join(test_path2, file)
file2 = 'sample0815.csv'
file_path3 = os.path.join(test_path2, file2)
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

    def test_trend_calc(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        key = ['SR520-Cal']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        d_test = dm.trend_calc(d_sample, key, min_size=5,
                               normalization='zscore',
                               method='pearsonr', visual=False)
        d_test1 = dm.trend_calc(d_sample, key, min_size=5,
                                normalization='log',
                                method='pearsonr', visual=False)
        d_test2 = dm.trend_calc(d_sample, key, min_size=5,
                                normalization='linear',
                                method='pearsonr', visual=False)
        assert len(d_test) > 0, 'no output'
        assert type(d_test) == pd.core.frame.DataFrame, 'wrong output'
        assert len(d_test1) > 0, 'no output'
        assert type(d_test1) == pd.core.frame.DataFrame, 'wrong output'
        assert len(d_test2) > 0, 'no output'
        assert type(d_test2) == pd.core.frame.DataFrame, 'wrong output'
        return

    def test_batch_alignment(self):
        d_model = pd.read_csv(file_path2)
        d_test = pd.read_csv(file_path3)
        d_merge = dm.batch_alignment(d_model, d_test)
        assert len(d_merge) > 0, 'alignment went wrong'
        assert type(d_merge) == pd.core.frame.DataFrame, 'wrong output'
        return

    def test_transpose(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        d_sample2 = dm.ms_cluster(d_sample, ['SR520-Cal'], 'linear',
                                  d_reduce=False, visual=False,
                                  eps=0.6, min_samples=10)
        d_ref = pd.read_csv(file_path3)
        d_model = d_sample2[d_sample2['label'] != -1]
        d_merge = dm.batch_alignment(d_model, d_ref)
        dilu_col = ['SR520-Cal', 'SR520_Cal']
        d_transpose = dm.transpose(d_merge, dilu_col)
        assert len(d_transpose) > 0, 'transpose went wrong'
        assert type(d_transpose) == pd.core.frame.DataFrame, 'wrong output'
        return

    def test_feature_model(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        d_sample2 = dm.ms_cluster(d_sample, ['SR520-Cal'], 'linear',
                                  d_reduce=False, visual=False,
                                  eps=0.6, min_samples=10)
        d_ref = pd.read_csv(file_path3)
        d_model = d_sample2[d_sample2['label'] != -1]
        d_merge = dm.batch_alignment(d_model, d_ref)
        dilu_col = ['SR520-Cal', 'SR520_Cal']
        d_transpose = dm.transpose(d_merge, dilu_col)
        reg = dm.feature_model(d_transpose, cluster_algorithm=False,
                               model_method='linear_reg', report=False)
        reg1 = dm.feature_model(d_transpose, cluster_algorithm=True,
                                model_method='ensemble_rf', report=False)
        reg2 = dm.feature_model(d_transpose, cluster_algorithm=False,
                                model_method='ensemble_bagging_dt',
                                report=False)
        reg3 = dm.feature_model(d_transpose, cluster_algorithm=True,
                                model_method='ensemble_bagging_svc',
                                report=False)
        assert type(reg) == sklearn.linear_model._base.LinearRegression, \
            'wrong regression applied'
        assert type(reg1) == list, 'wrong regression applied'
        assert type(reg2) == sklearn.ensemble._bagging.BaggingRegressor, \
            'wrong regression applied'
        assert type(reg3) == list, 'wrong regression applied'
        return

    def test_cluster_pred(self):
        keys = ['CEC', 'Blank', 'ISTD', 'Wash', 'Shutdown']
        d_ms = pd.read_csv(file_path)
        d_sample = dm.data_prep(d_ms, keys, rt_range=[1, 30],
                                mz_range=[200, 800], area_thres=500,
                                simp_summary=False)
        d_sample2 = dm.ms_cluster(d_sample, ['SR520-Cal'], 'linear',
                                  d_reduce=False, visual=False,
                                  eps=0.6, min_samples=10)
        d_ref = pd.read_csv(file_path3)
        d_model = d_sample2[d_sample2['label'] != -1]
        d_merge = dm.batch_alignment(d_model, d_ref)
        dilu_col = ['SR520-Cal', 'SR520_Cal']
        d_transpose = dm.transpose(d_merge, dilu_col)
        reg = dm.feature_model(d_transpose, cluster_algorithm=True,
                               model_method='ensemble_bagging_dt',
                               report=False)
        reg1 = dm.feature_model(d_transpose, cluster_algorithm=True,
                                model_method='ensemble_rf', report=False)
        dilu_test = dm.cluster_pred(reg, d_merge, '20181114_CoulterCreek_1',
                                    selected_cluster=[0, 1, 3, 5, 7, 9],
                                    model_merge='weighted_size')
        dilu_test2 = dm.cluster_pred(reg1, d_merge, '20181114_CoulterCreek_1',
                                     selected_cluster=[0, 1, 3, 5, 7, 9],
                                     model_merge='average')
        assert dilu_test >= 0, 'Invalid dilution prediction'
        assert dilu_test2 >= 0, 'Invalid dilution prediction'
        return
