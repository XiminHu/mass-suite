import pymzml
from tqdm import tqdm
import numpy as np

from scipy.integrate import simps
import pandas as pd

import peakutils

import glob
from pathlib import Path
import scipy
import pickle
import os
import re
import pyisopach
from scipy import special
import itertools
# from mss import mssdata
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
def noise_removal(mzml_scans, int_thres=1000):
    '''
    Remove mz&i pairs that i lower than int_thres
    from whole mzml file, looping through scans
    Only remove MS1 noise for now
    the output will overwrite the original mzml file
    int_thres: threshold for removing noises
    '''
    for scan in mzml_scans:
        if scan.ms_level == 1:
            drop_index = np.argwhere(scan.i <= int_thres)
            scan.i = np.delete(scan.i, drop_index)
            scan.mz = np.delete(scan.mz, drop_index)
        elif scan.ms_level != 1:
            continue

    return


# updated to all_than_close, when false only select closest one, when true
# append all, use as a backdoor for now if closest algorithm messed up
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


# def module?
this_dir, this_filename = os.path.split(__file__)
Model_file_t = os.path.join(this_dir, 'rfmodel_tuned.pkl')
rf_model_t = pickle.load(open(Model_file_t, 'rb'))
Pmodel = rf_model_t


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

        _, target_index = mz_locator(scan.mz, input_mz, error)
        if len(target_index) == 0:
            intensity.append(0)
        else:
            # intensity.append(scan.i[target_index])
            # CR -> if all_than_close=True
            # Change from sum to max
            intensity.append(max(scan.i[target_index]))

    return retention_time, intensity


def peak_pick(mzml_scans, input_mz, error, enable_score=True, peak_thres=0.01,
              peakutils_thres=0.02, min_d=1, rt_window=1.5,
              peak_area_thres=1e5, min_scan=5, max_scan=200, max_peak=5,
              overlap_tol=15, sn_detect=15):
    '''
    The function is used to detect peak for given m/z's chromatogram
    error: in ppm
    enable_score: option to enable the RF model
    peak_thres: base peak tolerance
    peakutils_thres: threshold from peakutils, may be repeated with peak_thres
    min_d: peaktuils parameter
    rt_window: window for integration only, didn't affect detection
    peak_area_thres: peak area limitation
    min_scan: min scan required to be detected as peak
    max_scan: max scan limit to exclude noise
    max_peak: max peak limit for selected precursor
    overlap_tot: overlap scans for two peaks within the same precursor
    sn_detect: scan numbers before/after the peak for sn calculation
    '''

    # Important funciont, may need to be extracted out later
    # Data output from the chromatogram_plot function

    rt, intensity = ms_chromatogram_list(mzml_scans, input_mz, error)

    # Get rt_window corresponding to scan number
    scan_window = int(
        (rt_window / (rt[int(len(intensity) / 2)] -
                      rt[int(len(intensity) / 2) - 1])) / 2)
    rt_conversion_coef = np.diff(rt).mean()
    # Get peak index
    indexes = peakutils.indexes(intensity, thres=peakutils_thres,
                                min_dist=min_d)

    result_dict = {}

    # dev note: boundary detection refinement
    for index in indexes:
        h_range = index
        l_range = index
        base_intensity = peak_thres * intensity[index]
        half_intensity = 0.5 * intensity[index]

        # Get the higher and lower boundary
        while intensity[h_range] >= base_intensity:
            h_range += 1
            if h_range >= len(intensity) - 1:
                break
            if intensity[h_range] < half_intensity:  # potentially record this
                if h_range - index > 4:  # fit r2 score,
                    # https://stackoverflow.com/questions/55649356/
                    # how-can-i-detect-if-trend-is-increasing-or-
                    # decreasing-in-time-series as alternative
                    x = np.linspace(h_range - 2, h_range, 3)
                    y = intensity[h_range - 2: h_range + 1]
                    (_slope, _intercept, r_value,
                     _p_value, _std_err) = scipy.stats.linregress(x, y)
                    # print(rt[h_range],r_value)
                    if abs(r_value) < 0.6:
                        break
        # Dev part 2, low priority since general peak shapes
        while intensity[l_range] >= base_intensity:
            l_range -= 1
            if l_range <= 1:
                break
            if intensity[l_range] < half_intensity:
                pass

        # Output a range for the peak list
        # If len(intensity) - h_range < 4:
        #     h_range = h_range + 3
        peak_range = []
        if h_range - l_range >= min_scan:
            if rt[h_range] - rt[l_range] <= rt_window:
                peak_range = intensity[l_range:h_range]
            else:
                if index - scan_window / 2 >= 1:
                    l_range = int(index - scan_window / 2)
                if index + scan_window / 2 <= len(intensity) - 1:
                    h_range = int(index + scan_window / 2)
                peak_range = intensity[l_range:h_range]
                # print(index + scan_window)

        # Follow Agilent S/N document
        width = rt[h_range] - rt[l_range]
        if len(peak_range) != 0:
            height = max(peak_range)
            hw_ratio = round(height / width, 0)
            neighbour_blank = (intensity[
                l_range - sn_detect: l_range] +
                intensity[h_range: h_range +
                          sn_detect + 1])
            noise = np.std(neighbour_blank)
            if noise != 0:
                sn = round(height / noise, 3)
            elif noise == 0:
                sn = 0

        # Additional global parameters
        # 1/2 peak range
        h_loc = index
        l_loc = index
        while intensity[h_loc] > half_intensity:
            h_loc += 1
            if h_loc >= len(intensity) - 1:
                break
        while intensity[l_loc] > half_intensity and l_loc > 0:
            l_loc -= 1

        # Intergration based on the simps function
        if len(peak_range) >= min_scan:
            integration_result = simps(peak_range)
            if integration_result >= peak_area_thres:
                # https://doi.org/10.1016/j.chroma.2010.02.010
                background_area = (h_range - l_range) * height
                ab_ratio = round(integration_result / background_area, 3)
                if enable_score is True:
                    h_half = h_loc + \
                        (half_intensity - intensity[h_loc]) / \
                        (intensity[h_loc - 1] - intensity[h_loc])
                    l_half = l_loc + \
                        (half_intensity - intensity[l_loc]) / \
                        (intensity[l_loc + 1] - intensity[l_loc])
                    # when transfer back use rt[index] instead
                    mb = (height - half_intensity) / \
                        ((h_half - index) * rt_conversion_coef)
                    ma = (height - half_intensity) / \
                        ((index - l_half) * rt_conversion_coef)
                    w = rt[h_range] - rt[l_range]
                    t_r = (h_half - l_half) * rt_conversion_coef
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

                # appending to result
                if len(result_dict) == 0:
                    (result_dict.update(
                     {index: [l_range, h_range,
                              integration_result, sn, score]}))
                # Compare with previous item
                elif integration_result != list(result_dict.values())[-1][2]:
                    s_window = abs(index - list(result_dict.keys())[-1])
                    if s_window > overlap_tol:
                        (result_dict.update(
                         {index: [l_range, h_range, integration_result,
                                  sn, score]}))
    # If still > max_peak then select top max_peak results
    if len(result_dict) > max_peak:
        result_dict = dict(sorted(result_dict.items(),
                                  key=lambda x: x[1][2], reverse=True))
        result_dict = dict(itertools.islice(result_dict.items(), max_peak))

    return result_dict


# Code review
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

    # Function to filter out empty mz slots to speed up the process
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

    mzlist = mz_gen(mzml_scans, err_ppm, mz_c_thres)
    print('Finding peaks...')

    result_dict = {}
    rt = []
    for scan in mzml_scans:
        rt.append(scan.scan_time[0])

    for mz in tqdm(mzlist):
        try:
            peak_dict = peak_pick(mzml_scans, mz, err_ppm, enable_score,
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


# Read in formula database
Formula_file = os.path.join(this_dir, '100-500.csv')
cfg = pd.read_csv(Formula_file, index_col=0)


def formula_calc(mz, composition, error=5, mode='pos'):
    '''
    Now only support 100-500 mz
    '''
    e_weight = 0.0005485799
    if mode == 'pos':
        mz = mz + e_weight
    elif mode == 'neg':
        mz = mz - e_weight

    low_mz = mz - mz * error * 1e-6
    high_mz = mz + mz * error * 1e-6

    hit = cfg[(cfg['Mass'] >= low_mz) & (cfg['Mass'] <= high_mz)].copy()
    hit['error in Da'] = mz - hit['Mass']
    hit['error in ppm'] = hit['error in Da'] / mz * 1e6

    hit = hit.iloc[hit['error in ppm'].abs().argsort()]

    for i in hit.index:
        f_reg = re.findall(r'([A-Z][a-z]*)', i)
        ip_reg = re.findall(r'([A-Z][a-z]*)', composition)
        if len([i for i in f_reg if i not in ip_reg]) > 0:
            hit.drop(i, inplace=True)

    return hit


def formula_prediction(mzml_scan, input_mz, error, composition='CHON',
                       f_error=5, mode='pos', relintensity_thres=1):
    '''
    Interactive spectrum plot with nearest retention time from the given scan
    mzml_scans: mzfile
    time: selected time for the scan
    '''

    def closest(lst, K):
        idx = np.abs(np.asarray(lst) - K).argmin()
        return idx

    intensity_max = ms_chromatogram_list(mzml_scan, input_mz, error)[1]
    scan = mzml_scan[np.argmax(intensity_max)]

    mz = scan.mz
    inten = scan.i

    precursor_idx = closest(mz, input_mz)
    precursor_mz = mz[precursor_idx]
    precursor_ints = inten[precursor_idx]

    rel_abundance = [i / precursor_ints * 100 for i in inten]

    prediction_table = formula_calc(
        precursor_mz, composition,
        error=f_error, mode='pos')

    # Find closest pair
    measured_spec = list(zip(mz, rel_abundance))

    def alpha(f):
        if f >= 80:
            alpha = 1
        elif 20 <= f < 80:
            alpha = -0.0033 * f + 1.2667
        elif 10 <= f < 20:
            alpha = -0.08 * f + 2.8
        elif 5 <= f < 10:
            alpha = -0.1 * f + 3
        elif 1 <= f < 5:
            alpha = -1.875 * f + 11.875
        return alpha

    def beta(f):
        if f >= 80:
            beta = 0.06
        elif 1 <= f < 80:
            beta = 2.0437 * (f ** 0.765)
        return beta

    for formula in prediction_table.index:
        mol = pyisopach.Molecule(formula)
        istp_mass = mol.isotopic_distribution()[0]
        istp_inten = mol.isotopic_distribution()[1]
        idx = np.argwhere(istp_inten >= relintensity_thres).reshape(1, -1)[0]
        # Ref SIRUIS 2013 paper
        m_list = istp_mass[idx]
        f_list = istp_inten[idx]

        theo_spec = list(zip(m_list, f_list))

        def dist(x, y):
            return (x[0]-y[0])**2 + (x[1]-y[1])**2

        score = []
        for p in theo_spec:

            measured_peak = [i for i in measured_spec if
                             i[0] >= p[0] * (1 - 1e-6 * 7)
                             and i[0] <= p[0] * (1 + 1e-6 * 7)]
            if len(measured_peak) != 0:
                hit_peak = min(measured_peak, key=lambda peak: dist(peak, p))
                # Function from SIRIUS, may need later validation
                f = hit_peak[1]
                sigma_m = np.sqrt(2) * hit_peak[0] * 1/3 * 1e-6 * 7 * alpha(f)
                x = abs(hit_peak[0] - p[0]) / sigma_m
                P_Mm = special.erfc(x)
                sigma_f = 1/3 * hit_peak[1] * np.log10(1 + beta(f))
                y = np.log10(hit_peak[1] / p[1]) / (np.sqrt(2) * sigma_f)
                P_fp = special.erfc(y)
                score.append(0.5 * (P_Mm + P_fp))
            else:
                hit_peak = []
                score.append(0)
        prediction_table.loc[formula, 'score'] = np.mean(score) * 100

    prediction_table.sort_values(by=['score'], inplace=True, ascending=False)

    return prediction_table
