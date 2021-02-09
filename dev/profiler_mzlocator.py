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
def mz_locator(input_list, mz, error, all_than_close=True):
    '''
    Find specific mzs from given mz and error range out from a given mz array
    input list: mz list
    mz: input_mz that want to be found
    error: error range is now changed to ppm level
    all_than_close: False only select closest one, True will append all
    '''
    target_mz = []
    target_index = []

    # ppm conversion
    error = error * 1e-6

    lower_mz = mz - error * mz
    higher_mz = mz + error * mz

    for i, mzs in enumerate(input_list):
        if mzs < lower_mz:
            continue
        elif mzs >= lower_mz:
            if mzs <= higher_mz:
                target_mz.append(mzs)
                target_index.append(i)

    if all_than_close is False:
        if len(target_mz) != 0:
            target_error = [abs(i - mz) for i in target_mz]
            minpos = target_error.index(min(target_error))
            t_mz = target_mz[minpos]
            t_i = target_index[minpos]
        else:
            t_mz = 0
            t_i = 'NA'
    if all_than_close is True:
        t_mz = target_mz
        t_i = target_index

    return t_mz, t_i
mz_locator(scans[1].mz,119.08,500)