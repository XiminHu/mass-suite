#!/usr/bin/env python3
# test.py
import time
import sys
# import mss
sys.path.append('../')
from mss import visreader as mvis
from mss import mssmain as msm
from mss import align
import pandas as pd
import numpy as np
from tqdm import tqdm
import peakutils
import scipy
from scipy.integrate import simps
import itertools
path = '../example_data/ex_1.mzML'
scans = msm.get_scans(path, ms_all=False, ms_lv=1)
#noise removal
msm.noise_removal(scans, 2000)
scans = scans[:200]

@profile
def mz_gen(mzml_scans, err_ppm, mz_c_thres):
    # Function remake needed
    pmz = []
    for scan in mzml_scans:
        pmz.append(scan.mz)
    pmz = np.hstack(pmz).squeeze()

    # According to msdial it should be mz + error * mz
    # To avoid mz slicing issue
    # Gap used to be 2*error*mz
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4449330/#SD1
    def mz_list_gen(minmz, maxmz, error_ppm):
        error = error_ppm * 1e-6
        mz_list = [minmz]
        mz = minmz
        while mz <= maxmz:
            mz = mz + error * mz
            mz_list.append(mz)
        return mz_list

    mz_list = mz_list_gen(pmz.min(), pmz.max(), err_ppm)

    final_mz = []
    for m in mz_list:
        lm = m - err_ppm * 1e-6 * m
        hm = m + err_ppm * 1e-6 * m
        if len(pmz[(pmz <= hm) & (pmz >= lm)]) >= mz_c_thres:
            final_mz.append(m)

    return final_mz

mz_gen(scans,20,5)