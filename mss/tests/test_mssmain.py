import os
import unittest
import numpy as np
import mss
from mss import mssmain

data_path = os.path.join(mss.__path__[0], 'data')


class test_mssmain(unittest.TestCase):

    def test_get_scans(self):
        """Tests stack functionlity"""
        path = '../../example_data/ex_1.mzML'
        test_scan = mssmain.get_scans(path)
        assert type(test_scan) == list, "output type is wrong"
        assert len(test_scan) == 4061, "output length is wrong"
        assert len(test_scan[0].mz) != 0, "mz array didn't read in"
        assert len(test_scan[0].i) != 0, "i array didn't read in"
        return

    def test_noise_removal(self):
        path = '../../example_data/ex_1.mzML'
        test_scan = mssmain.get_scans(path)
        test_scan = mssmain.noise_removal(test_scan)
        assert np.where(test_scan[0].i < 5000) == 0,\
            "noise didn't properly removed"
        return

    def test_mzlocator(self):
        path = '../../example_data/ex_1.mzML'
        test_scan = mssmain.get_scans(path)
        mzlist = mssmain.mz_locator(test_scan[0].mz, 117.113, 20)
        assert type(mzlist) == tuple, "output type is wrong"
        assert mzlist[0][0] - 117.113 <= 20 * 117.113,\
            "selected mz is above defined error range"
        assert mzlist[1][0] != 0, "index error"
        return

    def test_peak_pick(self):  # More to add during later development
        path = '../../example_data/ex_1.mzML'
        test_scan = mssmain.get_scans(path)
        test_dict = mssmain.peak_pick(test_scan, 299.146, 20)
        assert type(test_dict) == dict, "output type is wrong"
        assert len(test_dict.keys()) != 0, "didn't find the peak"
        assert list(test_dict.values())[0][2] - 136350 <= 1,\
            "intergration result error"
        return

    def test_peak_list(self):  # More to add during later development
        path = '../../example_data/ex_1.mzML'
        test_scan = mssmain.get_scans(path)
        d_test = mss.peak_list(test_scan[:200], 20)
        assert d_test.shape == (47, 5), "feature seeking is off"
        return

    def test_batch_scans(self):
        return

    def test_batch_peak(self):
        return
