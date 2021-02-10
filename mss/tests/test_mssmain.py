import os
import unittest
import mss
from mss import mssmain

test_path = os.path.join(mss.__path__[0], 'tests')
data_path = os.path.join(test_path, 'data')
file = 'ex_1.mzML'
file_path = os.path.join(data_path, file)
test_scan_file = mssmain.get_scans(file_path)


class test_mssmain(unittest.TestCase):

    def test_get_scans(self):
        """Tests stack functionlity"""
        test_scan = mssmain.get_scans(file_path)
        assert type(test_scan) == list, "output type is wrong"
        assert len(test_scan) == 4061, "output length is wrong"
        assert len(test_scan[0].mz) != 0, "mz array didn't read in"
        assert len(test_scan[0].i) != 0, "i array didn't read in"
        return

    def test_noise_removal(self):
        test_scan = test_scan_file
        mssmain.noise_removal(test_scan)
        assert sum(test_scan[0].i <= 1000) == 0,\
            "noise didn't properly removed"
        return

    def test_mzlocator(self):
        test_scan = test_scan_file
        mzlist = mssmain.mz_locator(test_scan[0].mz, 117.113, 20)
        assert type(mzlist) == tuple, "output type is wrong"
        assert mzlist[0][0] - 117.113 <= 20 * 117.113,\
            "selected mz is above defined error range"
        assert mzlist[1][0] != 0, "index error"
        return

    def test_peak_pick(self):  # More to add during later development
        test_scan = test_scan_file
        test_dict = mssmain.peak_pick(test_scan, 299.146, 20)
        assert type(test_dict) == dict, "output type is wrong"
        assert len(test_dict.keys()) != 0, "didn't find the peak"
        assert list(test_dict.values())[0][2] - 136350 <= 1,\
            "intergration result error"
        return

    def test_peak_list(self):  # More to add during later development
        test_scan = test_scan_file
        d_test = mssmain.peak_list(test_scan[:200], 20, enable_score=True)
        assert d_test.shape[1] == 5, "score column is off"
        assert d_test.score.dtype == 'int64', "wrong score type"
        return

    def test_formula_calc(self):  # More to add during later development
        test_mz = 299.1765
        test_composition = 'CHON'
        d_test = mssmain.formula_calc(test_mz, test_composition)
        assert d_test.shape[0] == 2, "Wrong shape"
        assert d_test.index.dtype == 'object', "Wrong dtype"
        return

    def test_formula_prediction(self):
        test_scan = test_scan_file
        d_test = mssmain.formula_prediction(test_scan, 299.1765, 5, f_error=10)
        assert d_test.shape[1] == 4, "Wrong dataframe shape"
        return


#    def test_batch_scans(self):
#       return

#    def test_batch_peak(self):
#        return
