import pymzml
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np

import pyteomics
from pyteomics import mzml, auxiliary
import matplotlib.pyplot as plt
import numpy as np
import math
import plotly.graph_objects as go
import re
from scipy.integrate import simps
import pandas as pd

import peakutils
from peakutils.plot import plot as pplot


def get_scans(path, ms_lv = 1):
    run = pymzml.run.Reader(path)
    scans = []
    for scan in run:
        if scan.ms_level == ms_lv:
            scans.append(scan)
    return scans


#Noise removal
def noise_removal(mzml_scans, int_thres = 5000):
    '''
    the output will overwrite the original mzml file
    '''
    for scan in mzml_scans:
        drop_index = np.argwhere(scan.i <= int_thres)
        scan.i = np.delete(scan.i, drop_index)
        scan.mz = np.delete(scan.mz, drop_index)
    
    return



def mz_locator(input_list, mz, error):
    '''
    Find specific mzs from given mz and error range
    input list: mz list
    '''
    target_mz = []
    target_index = []
    
    lower_mz = mz - error
    higher_mz = mz + error

    for i, mzs in enumerate(input_list):
        if mzs < lower_mz:
            continue
        elif mzs >= lower_mz:
            if mzs <= higher_mz:
                target_mz.append(mzs)
                target_index.append(i)
        elif mzs > higher_mz:
                target_mz = 0
                target_index = 'NA'
                break
        
    return target_mz, target_index


#Deal with roi, check peakonly
def peak_pick(mzml_scans, input_mz, error, peak_base = 5000, thr = 0.02, min_d = 1, rt_window = 1.5, peak_area_thres = 1e5, min_scan = 15, max_scan = 200, max_peak = 7):
    '''
    rt, ints from ms_chromatogram_list
    rt_window now set up for minutes
    '''
    
    #Important funciont, may need to be extracted out later
    def ms_chromatogram_list(mzml_scans, input_mz, error):
        '''
        Generate a peak list for specific input_mz over whole rt period from the mzml file
        ***Most useful function!
        '''
        retention_time = []
        intensity = []
        for scan in mzml_scans:
            #print(i)
            retention_time.append(scan.scan_time[0])

            target_mz, target_index = mz_locator(scan.mz, input_mz, error)
            if target_index == 'NA':
                intensity.append(0)
            else:
                intensity.append(sum(scan.i[target_index]))
            
        return retention_time, intensity
    
    rt, intensity = ms_chromatogram_list(mzml_scans, input_mz, error)
    
    #Get rt_window corresponded scan number
    scan_window = int((rt_window / (rt[int(len(intensity) / 2)] - rt[int(len(intensity) / 2) - 1])) / 2)
    
    #Get peak index
    indexes = peakutils.indexes(intensity, thres=thr, min_dist = min_d)
    
    result_dict = {}
    
    for index in indexes:
        h_range = index
        l_range = index
        base_intensity = peak_base

        #Get the higher and lower boundary
        while intensity[h_range] >= base_intensity:
            h_range += 1
            if h_range > len(intensity)-2:
                break
        while intensity[l_range] >= base_intensity:
            l_range -= 1
        #Output a range from the peak list
        
        peak_range = []
        if h_range - l_range >= min_scan:
            if rt[h_range] - rt[l_range] <= rt_window:
                peak_range = intensity[l_range:h_range]
            else:
                l_range = index - scan_window
                h_range = index + scan_window
                peak_range = intensity[l_range:h_range]
                #print(index + scan_window)

        #Intergration based on the simps function
        if len(peak_range) >= min_scan:
            integration_result = simps(peak_range)
            if integration_result >= peak_area_thres:
                if len(result_dict) == 0:
                    result_dict.update({index : [l_range, h_range, integration_result]})
                elif integration_result != list(result_dict.values())[-1][2]: #Compare with previous item
                    result_dict.update({index : [l_range, h_range, integration_result]})
                
        #Filtering:
        #1. delete results that l_range/h_range within 5 scans
        #3. If still >5 then select top 5 results
        #list(result_dict.values())[-1]
    
    #Noise filter
    if len(result_dict) > max_peak:
        result_dict = {}
        


    return result_dict


#Main dev issue--iteration efficiency & peak picking algorithm update
#2way: 1. optimize iteration algorithm, 2. optimize mz_list selection--drop useless mzs
def peak_list(mzml_scans, mz_error, mzmin = 0, mzmax=1000, peak_base = 0.005, thr = 0.02, min_d = 1, rt_window = 1.5, peak_area_thres = 1e5, min_scan = 7, scan_thres = 7):
    '''
    input from ms1_baseline function
    Q to solve: how to correctly select mz slice?? see mz_locator
    '''
    
    #Get m/z range
    def mz_array(scans, mz_error):
        min_mz = scans[0].mz.min()
        max_mz = scans[0].mz.max()
        for scan in scans:
            if min_mz > scan.mz.min():
                min_mz = scan.mz.min()
            if max_mz < scan.mz.max():
                max_mz = scan.mz.max()

        mz_list = np.arange(min_mz + mz_error, max_mz - mz_error,2 * mz_error)

        return mz_list
    
    mzlist = mz_array(mzml_scans, mz_error)
    
    if mzlist.max() >= mzmax:
        mzlist = mzlist[mzlist <= mzmax]
    if mzlist.min() <= mzmin:
        mzlist = mzlist[mzlist >= mzmin]
    
    
    result_dict = {}
    rt = []
    for scan in mzml_scans:
        rt.append(scan.scan_time[0])
    
    for mz in tqdm(mzlist):
        try:
            peak_dict = peak_pick(mzml_scans, mz, mz_error)
        except:
            pass
        
        if len(result_dict) == 0:
            for index in peak_dict:
                result_dict.update({'m/z' : [mz],
                                    'rt' : [rt[index]],
                                    'peak area' : [peak_dict[index][2]]})
        else:
            for index in peak_dict:
                result_dict['m/z'].append(mz)
                result_dict['rt'].append(rt[index])
                result_dict['peak area'].append(peak_dict[index][2])
            
    print('Peak processing finished!')
    d_result = pd.DataFrame(result_dict)
    print('Dataframe created!')
    
    return d_result