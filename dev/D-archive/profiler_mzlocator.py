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

@profile
def boolean_index(input_list, mz, error):
    array = np.asarray(input_list)
    error = error * 1e-6

    lower_mz = mz - error * mz
    higher_mz = mz + error * mz
    

    index = (array >= lower_mz) & (array <= higher_mz)

    return array[index],np.where(index)[0]

print('Boolean index:\t', end='')

boolean_index(scans[1].mz,119.08,500)