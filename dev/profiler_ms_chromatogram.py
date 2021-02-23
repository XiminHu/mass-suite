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
def ms_chromatogram_list(mzml_scans, input_mz, error):
    '''
    Generate a peak list for specific input_mz over
    whole rt period from the mzml file
    ***Most useful function!
    '''
    intensity = []
    retention_time = [i.scan_time[0] for i in mzml_scans]
    for scan in mzml_scans:
        _, target_index = msm.mz_locator(scan.mz, input_mz, error)
        if len(target_index) == 0:
            intensity.append(0)
        else:
            # intensity.append(scan.i[target_index])
            # CR -> if all_than_close=True
            # Change from sum to max
            intensity.append(max(scan.i[target_index]))

    return retention_time, intensity

ms_chromatogram_list(scans,299.1765,500)