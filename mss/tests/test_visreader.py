import unittest
from mss import visreader


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
