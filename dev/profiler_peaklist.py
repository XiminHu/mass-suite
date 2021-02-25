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
path = 'D:/UW/mssvalidation/20201119_ENTACT_validation_data/20201117_ENTACT_505_1.mzML'
scans = msm.get_scans(path, ms_all=False, ms_lv=1)
#noise removal
msm.noise_removal(scans, 2000)
scans1 = scans[:100]

@profile
def peak_list(mzml_scans, err_ppm=10, enable_score=True, mz_c_thres=5,
              peak_base=0.005, peakutils_thres=0.02, min_d=1, rt_window=1.5,
              peak_area_thres=1e5, min_scan=5, max_scan=50,
              max_peak=5):
    '''
    Generate a dataframe by looping throughout the
    whole mz space from a given mzml file
    ref to peak_picking function
    all the parameters included in peak_pick
    mz_c_thres: defines how much mz need to be within a cluster for
    a valid precursor in peak list detection
    '''

    # Get m/z range -- updated 0416
    print('Generating mz list...')

    mzlist = msm.mz_gen(mzml_scans, err_ppm, mz_c_thres)
    print('Finding peaks...')

    result_dict = {}
    rt = []
    for scan in mzml_scans:
        rt.append(scan.scan_time[0])

    for mz in tqdm(mzlist):
        # * python instrumentation run time
        # * cython to rewrite
        try:
            peak_dict = msm.peak_pick(mzml_scans, mz, err_ppm, enable_score,
                                  peak_thres=peak_base,
                                  peakutils_thres=peakutils_thres,
                                  min_d=min_d, rt_window=rt_window,
                                  peak_area_thres=peak_area_thres,
                                  min_scan=min_scan, max_scan=max_scan,
                                  max_peak=max_peak)
        except Exception:  # Catch exception?
            peak_dict = {}

        if len(peak_dict) != 0:
            if len(result_dict) == 0:
                for index in peak_dict:
                    result_dict.update({'m/z': [mz],
                                        'rt': [rt[index]],
                                        'sn': [peak_dict[index][3]],
                                        'score': [peak_dict[index][4]],
                                        'peak area': [peak_dict[index][2]]})
            else:
                for index in peak_dict:
                    result_dict['m/z'].append(mz)
                    result_dict['rt'].append(rt[index])
                    result_dict['sn'].append(peak_dict[index][3])
                    result_dict['score'].append(peak_dict[index][4])
                    result_dict['peak area'].append(peak_dict[index][2])
    # print(result_dict)
    print('Peak processing finished!')
    d_result = pd.DataFrame(result_dict)
    d_result['rt'] = round(d_result['rt'], 2)
    d_result['m/z'] = round(d_result['m/z'], 4)
    print('Dataframe created!')

    return d_result

peak_list(scans[:500], err_ppm=10,enable_score=False)