import unittest
from mss import visreader
import os
import mss
from mss import mssmain
data_path = os.path.dirname(os.path.join(mss.__path__[0]))
batchpath = data_path + '/example_data/ex_1.mzML'


class test_dm(unittest.TestCase):
    def test_mz_locator(self):
        test_array = [0, 10, 20, 50, 100]
        result = visreader.mz_locator(test_array, 49.9999, 20)
        assert result[0][0] == 50, 'Wrong mz found'
        assert result[1][0] == 3, 'Wrong index found'
        return

    def test_formula_mass(self):
        test_formula = 'C18H22N2O2'
        result = visreader.formula_mass(test_formula)
        assert result == 299.1764085799, 'Wrong output'
        return

    def test_tic_plot(self):
        scans = mssmain.get_scans(batchpath, ms_all=False, ms_lv=1)
        test_plt = visreader.tic_plot(scans, interactive=False)
        assert isinstance(test_plt, type(None)), 'invalid plot'
        return

    def test_ms_plot(self):
        scans = mssmain.get_scans(batchpath, ms_all=False, ms_lv=1)
        test_plt = visreader.ms_plot(scans, 8.47, True)
        assert isinstance(test_plt, type(None)), 'invalid result'
        return

    def test_ms_chromatogram(self):
        scans = mssmain.get_scans(batchpath, ms_all=False, ms_lv=1)
        test_plt = visreader.ms_chromatogram(scans, 'C18H22N2O2',
                                             20, False, 'pos',
                                             False, False, 'pubchem')
        assert isinstance(test_plt, type(None)), 'invalid result'
        return
