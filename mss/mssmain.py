import pymzml
from tqdm import tqdm
import numpy as np

from scipy.integrate import simps
import pandas as pd

import peakutils

import glob
from pathlib import Path
import scipy
from mss import mssdata
# Modeling modules
# from tensorflow import keras
# import h5py
# from mssdata import peakmodel


def get_scans(path, ms_all=False, ms_lv=1):
    '''
    The function is used to reorganize the pymzml reading
    into a list that will have better access
    path: input mzml path
    ms_all: if you want all the ms_level be check
    ms_lv: ms_level you want to export
    '''

    # Read path using pymzml
    run = pymzml.run.Reader(path)

    scans = []
    if ms_all is False:
        for scan in run:
            if scan.ms_level == ms_lv:
                scans.append(scan)
    elif ms_all is True:
        for scan in run:
            scans.append(scan)

    return scans


# Noise removal
def noise_removal(mzml_scans, int_thres=5000):
    '''
    Remove mz&i pairs that i lower than int_thres
    from whole mzml file, looping through scans
    the output will overwrite the original mzml file
    '''
    for scan in mzml_scans:
        drop_index = np.argwhere(scan.i <= int_thres)
        scan.i = np.delete(scan.i, drop_index)
        scan.mz = np.delete(scan.mz, drop_index)

    return


# Code review
# updated to select_app, when false only select closest one, when true
# append all, use as a backdoor for now if closest algorithm messed up
def mz_locator(input_list, mz, error, select_app=True):
    '''
    Find specific mzs from given mz and error range out from a given mz array
    input list: mz list
    mz: input_mz that want to be found
    error: error range is now changed to ppm level
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

    if select_app is False:
        if len(target_mz) != 0:
            target_error = [abs(i - mz) for i in target_mz]
            minpos = target_error.index(min(target_error))
            t_mz = target_mz[minpos]
            t_i = target_index[minpos]
        else:
            t_mz = 0
            t_i = 'NA'
    if select_app is True:
        t_mz = target_mz
        t_i = target_index

    return t_mz, t_i

# Read model for peak assessment


Pmodel = mssdata.peakmodel.rf_model_t


def peak_pick(mzml_scans, input_mz, error, enable_score=True, peak_thres=0.01,
              thr=0.02, min_d=1, rt_window=1.5, peak_area_thres=1e5,
              min_scan=15, max_scan=200, max_peak=5,
              min_scan_window=20, sn_range=7):
    '''
    firstly get rt, intensity from given mz and error out of the mzml file
    Then find peak on the intensity array, represent as index --> index
    Find peak range by looping from peak index forward/backward
    until hit the peak_base --> l_range,h_range. peakspan = h_range - l_range
    Trim/correct peak range is too small or too large,
    using min_scan/max_scan,min_scan_window --> trimed l/h_range
    Integration of peak based on the given range
    using simp function --> peakarea
    '''

    # Important funciont, may need to be extracted out later
    # Data output from the chromatogram_plot function
    def ms_chromatogram_list(mzml_scans, input_mz, error):
        '''
        Generate a peak list for specific input_mz over
        whole rt period from the mzml file
        ***Most useful function!
        '''

        # Create empty list to store the data
        retention_time = []
        intensity = []
        for scan in mzml_scans:
            # print(i)
            retention_time.append(scan.scan_time[0])

            target_mz, target_index = mz_locator(scan.mz, input_mz, error)
            if len(target_index) == 0:
                intensity.append(0)
            else:
                # intensity.append(scan.i[target_index])
                # if select_app=True
                intensity.append(sum(scan.i[target_index]))

        return retention_time, intensity

    rt, intensity = ms_chromatogram_list(mzml_scans, input_mz, error)

    # Get rt_window corresponded scan number
    scan_window = int(
                      (rt_window / (rt[int(len(intensity) / 2)] -
                       rt[int(len(intensity) / 2) - 1])) / 2)
    rt_conversion_rate = rt[1] - rt[0]
    # Get peak index
    indexes = peakutils.indexes(intensity, thres=thr, min_dist=min_d)

    result_dict = {}

    # dev note: boundary detection refinement
    for index in indexes:
        h_range = index
        l_range = index
        # use relative thres, also considering S/N, 1/2 rt point?
        base_intensity = peak_thres * intensity[index]
        half_intensity = 0.5 * intensity[index]

        # Get the higher and lower boundary
        while intensity[h_range] >= base_intensity:
            h_range += 1
            if intensity[h_range] < half_intensity:  # potentially record this
                if h_range - index > 4:  # fit r2 score,
                    # keep record
                    # https://stackoverflow.com/questions/55649356/
                    # how-can-i-detect-if-trend-is-increasing-or-
                    # decreasing-in-time-series as alternative
                    x = np.linspace(h_range - 2, h_range, 3)
                    y = intensity[h_range - 2: h_range + 1]
                    (slope, intercept, r_value,
                     p_value, std_err) = scipy.stats.linregress(x, y)
                    # print(rt[h_range],r_value)
                    if abs(r_value) < 0.6:
                        break
                    elif h_range > len(intensity) - 2:
                        break
        # Dev part 2, low priority since general peak shapes
        while intensity[l_range] >= base_intensity:
            l_range -= 1
            if intensity[l_range] < half_intensity:
                pass  # backdoor for recording 1/2 rt point
        # Output a range from the peak list

        peak_range = []
        if h_range - l_range >= min_scan:
            if rt[h_range] - rt[l_range] <= rt_window:
                peak_range = intensity[l_range:h_range]
            else:
                l_range = index - scan_window
                h_range = index + scan_window
                peak_range = intensity[l_range:h_range]
                # print(index + scan_window)

        # Calculate for S/N
        signal = intensity[index]
        neighbour_blank = (intensity[
                           l_range - sn_range: l_range] +
                           intensity[h_range + 1: h_range +
                           sn_range + 1])
        noise = max(neighbour_blank)
        if noise != 0:
            sn = round(signal / noise, 3)
        elif noise == 0:
            sn = 0

        # Calculate height/width, consider log10 transform
        height = signal
        width = rt[h_range] - rt[l_range]
        hw_ratio = round(height / width, 0)

        # Additional global parameters
        # 1/2 peak range
        h_loc = index
        l_loc = index
        while intensity[h_loc] > half_intensity:
            h_loc += 1
        while intensity[l_loc] > half_intensity:
            l_loc -= 1
        # calculate for slope -- interpolation included-- pay attention!
        h_half = h_loc + \
            (half_intensity - intensity[h_loc]) / \
            (intensity[h_loc - 1] - intensity[h_loc])
        l_half = l_loc + \
            (half_intensity - intensity[l_loc]) / \
            (intensity[l_loc + 1] - intensity[l_loc])
        # when transfer back use rt[index] instead
        mb = (height - half_intensity) / \
            ((h_half - index) * rt_conversion_rate)
        ma = (height - half_intensity) / \
            ((index - l_half) * rt_conversion_rate)

        # Intergration based on the simps function
        if len(peak_range) >= min_scan:
            integration_result = simps(peak_range)
            if integration_result >= peak_area_thres:
                # Calculate Area/background ratio, i.e, peak area vs
                # rectangular area as whole(if =1 then peak is a pleateu)
                background_area = (h_range - l_range) * height
                ab_ratio = round(integration_result / background_area, 3)

                # Awaiting to be added:
                # model prediction as the assessment score column
                # score = model.fit(X_para from all above)
                if enable_score is True:
                    w = rt[h_range] - rt[l_range]
                    t_r = (h_half - l_half) * rt_conversion_rate
                    l_width = rt[index] - rt[l_range]
                    r_width = rt[h_range] - rt[index]
                    assym = r_width / l_width
                    var = (w ** 2 / (1.764 * ((r_width / l_width)
                           ** 2) - 11.15 * (r_width / l_width) + 28))
                    x_peak = [w, t_r, l_width, r_width, assym,
                              integration_result, sn, hw_ratio, ab_ratio,
                              height, ma, mb, ma + mb, mb / ma, var]
                    x_input = np.asarray(x_peak)
                    # score = np.argmax(Pmodel.predict(x_input.reshape(1,-1)))
                    # for tensorflow
                    score = int(Pmodel.predict(x_input.reshape(1, -1)))
                elif enable_score is False:
                    score = 1
                # final result append score

                # appending to result
                if len(result_dict) == 0:
                    (result_dict.update(
                     {index: [l_range, h_range,
                      integration_result, sn, score]}))
                # Compare with previous item
                elif integration_result != list(result_dict.values())[-1][2]:
                    s_window = abs(index - list(result_dict.keys())[-1])
                    if s_window > min_scan_window:
                        (result_dict.update(
                         {index: [l_range, h_range, integration_result,
                          sn, score]}))

        # Filtering:
        # 1. delete results that l_range/h_range within 5 scans
        # 3. If still >5 then select top 5 results
        # list(result_dict.values())[-1]

    # Noise filter
    if len(result_dict) > max_peak:
        result_dict = {}

    return result_dict


# Code review
def peak_list(mzml_scans, err_ppm=20, enable_score=True, mz_c_thres=5,
              peak_base=0.005, thr=0.02, min_d=1, rt_window=1.5,
              peak_area_thres=1e5, min_scan=7, scan_thres=7):
    '''
    Generate a dataframe by looping throughout the
    whole mz space from a given mzml file
    ref to peak_picking function
    Q to solve: how to correctly select mz slice?? see mz_locator
    '''

    # Get m/z range -- updated 0416
    print('Generating mz list...')

    # Function to filter out empty mz slots to speed up the process
    def mz_gen(mzml_scans, err_ppm, mz_c_thres):
        pmz = []
        for scan in mzml_scans:
            pmz.append(scan.mz)
        pmz = np.hstack(pmz).squeeze()

        # Function to generate a reference mz list using a defined step
        # according to user setting
        def mz_list_gen(minmz, maxmz, error_ppm):
            error = error_ppm * 1e-6
            mz_list = [minmz]
            mz = minmz
            while mz <= maxmz:
                mz = mz + 2 * error * mz
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

    mzlist = mz_gen(mzml_scans, err_ppm, mz_c_thres)  # New list for looping
    print('Finding peaks...')

    result_dict = {}
    rt = []
    for scan in mzml_scans:
        rt.append(scan.scan_time[0])

    for mz in tqdm(mzlist):
        try:
            peak_dict = peak_pick(mzml_scans, mz, err_ppm, enable_score)
        except Exception:
            pass

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


# Only work on MS1 scans, needs update on the MS2 included scans
def batch_scans(path, remove_noise=True, thres_noise=1000):
    all_files = glob.glob(path + "/*.mzML")
    scans = []
    file_list = []
    for file in tqdm(all_files):
        scan = get_scans(file)
        if remove_noise is True:
            noise_removal(scan, thres_noise)
        scans.append(scan)
        file_list.append(Path(file).name)
    print(file_list)
    print('Batch read finished!')

    return scans, file_list


def batch_peak(batch_input, source_list, mz, error):
    rt_max = []
    rt_start = []
    rt_end = []
    peak_area = []
    source = []
    for i, scans in enumerate(batch_input):
        rt = []
        result_dict = peak_pick(scans, mz, error)
        for scan in scans:
            rt.append(scan.scan_time[0])
        for index in result_dict:
            rt_max.append(round(rt[index], 2))
            rt_start.append(round(rt[list(result_dict.values())[0][0]], 2))
            rt_end.append(round(rt[list(result_dict.values())[0][1]], 2))
            peak_area.append(round(result_dict[index][2], 4))
            source.append(source_list[i])
    result_dict = {'rt_max': rt_max,
                   'rt_start': rt_start,
                   'rt_end': rt_end,
                   'peak_area': peak_area,
                   'source': source
                   }
    d_result = pd.DataFrame(result_dict)

    return d_result
