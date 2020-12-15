import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from mssmain import ms_chromatogram_list, batch_scans, peak_list
"""This file contains the functions needed to correct the peak positions
of different samples accumulated in one day's worth of analysis."""


def stack(d_batch):
    """This function compiles all samples files into one dataframe
    for analysis. It takes in files that are of the .txt type."""
    all_samples = []
    # Calculating the number of files in path
    num_files = len(d_batch)
    for i in range(num_files):
        sample_name = d_batch[i]
        sample_df = pd.DataFrame(sample_name, columns=['m/z',
                                 'rt', 'sn', 'score', 'peak area'],
                                 dtype=np.float32)
        all_samples += [sample_df]
    # Combining all the dataframes into one
    total_samples = pd.concat(all_samples)
    all_samples.clear()
    # Cleaning up the data before processing
    total_samples.loc[total_samples.sn == 0, 'sn'] = float('inf')
    total_samples.loc[total_samples.score == 3.0, 'score'] = 0.1
    total_samples.loc[total_samples.score == 2.0, 'score'] = 0.6
    print("Process completed!")
    return num_files, total_samples


def realignment(d_batch, export_name, name_list, align_rt_err, align_mz_err):
    """This function works by using one .txt file as a reference in which
    other files realigned to in terms of precursor and RT. """
    np.warnings.filterwarnings('ignore') # ignore warning from sn inf value
    RT_error = align_rt_err  # units of minutes, can be adjusted
    alignment_df = pd.DataFrame()
    mz_error = align_mz_err
    # standard_sample = os.listdir(batch_path)[0]  # first sample_name
    # reads .txt file into a dataframe
    # standard_df = pd.read_csv(batch_path + standard_sample, usecols=['rt',
    #                          'm/z', 'sn', 'score', 'peak area'],
    #                          dtype=np.float32)
    standard_df = pd.DataFrame(d_batch[0], columns=['m/z',
                            'rt', 'sn', 'score', 'peak area'],
                            dtype=np.float32)
    # Creating the reference df based off 1st sample_name
    alignment_df['m/z'] = standard_df.iloc[:, 0]
    alignment_df['rt'] = standard_df.iloc[:, 1]
    alignment_df['sn'] = standard_df.iloc[:, 2]
    alignment_df['score'] = standard_df.iloc[:, 3]
    alignment_df['Sum RT (min)'] = 0.0
    alignment_df['Sum Precursor m/z'] = 0.0
    alignment_df['Sum sn'] = 0.0
    alignment_df['Sum score'] = 0.0
    alignment_df['Sum RT (min)'] = alignment_df[
                                   'Sum RT (min)'].astype('float32')
    alignment_df['Sum Precursor m/z'] = alignment_df[
                                        'Sum Precursor m/z'].astype('float32')
    alignment_df['Sum sn'] = alignment_df['Sum sn'].astype('float32')
    alignment_df['Sum score'] = alignment_df['Sum score'].astype('float32')
    alignment_df.loc[alignment_df.sn == 0, 'sn'] = float('inf')
    alignment_df.loc[alignment_df.score == 3.0, 'score'] = 0.1
    alignment_df.loc[alignment_df.score == 2.0, 'score'] = 0.6
    # col_index to track where to add sample_name columns
    col_index = 8
    num_files, total_samp = stack(d_batch)
    for files in range(num_files):
        sample_name = name_list[files]  # Creates file columns
        sample_name = sample_name[:-5]  # removes extension part in sample_name name
        alignment_df[sample_name] = 0.0
        alignment_df[sample_name] = alignment_df[sample_name].astype('float32')
    print("Initial reference built")
    print("Alignment beginning..")
    for row in tqdm(range(len(total_samp))):
        if alignment_df.isnull().values.any() is True:
            alignment_df = alignment_df.fillna(0)
        else:
            pass
        row_max = len(alignment_df)
        if row < range(len(total_samp))[-1]:
            # Loops until a change in sample_name occurs by index
            if total_samp.index[row] < total_samp.index[row + 1]:
                # checks to see if row exists in reference or not
                overlap = np.where(((alignment_df.iloc[:, 1] - RT_error)
                                   <= total_samp.iloc[row, 1]) &
                                   (total_samp.iloc[row, 1] <=
                                   (alignment_df.iloc[:, 1] + RT_error))
                                   & ((alignment_df.iloc[:, 0] - mz_error)
                                      <= total_samp.iloc[row, 0])
                                   & (total_samp.iloc[row, 0] <=
                                       (alignment_df.iloc[:, 0] + mz_error)))
                if len(overlap[0]) > 0:
                    overlap = [overlap][0][0][0]
                    # checks for any duplicates
                    if (alignment_df.iloc[overlap, col_index]) > 0:
                        alignment_df.at[overlap, alignment_df.columns[
                                        col_index]] += total_samp.iloc[row, 4]
                        # Averages the mean if there is a repeat intensity
                        alignment_df.at[overlap, alignment_df.columns[
                                         col_index]] /= 2
                    else:
                        # sets intensity value where realignment was identified
                        alignment_df.at[overlap, alignment_df.columns[
                                             col_index]] = total_samp.iloc[
                                                           row, 4]
                        # adds rt, sn, m/z, score to aligned point
                        alignment_df.at[overlap,
                                        'Sum RT (min)'] += (total_samp.
                                                            iloc[row, 1])
                        alignment_df.at[overlap,
                                        'Sum Precursor m/z'] += (total_samp.
                                                                 iloc[row, 0])
                        alignment_df.at[overlap, 'Sum sn'] += (total_samp.
                                                               iloc[row, 2])
                        alignment_df.at[overlap, 'Sum score'] += (total_samp.
                                                                  iloc[row, 3])
                else:
                    # appends as a new reference if nothing is found
                    alignment_df.loc[row_max] = total_samp.iloc[row, 0:4]
                    alignment_df.at[alignment_df.index[
                                   -1], alignment_df.columns[
                                        col_index]] = total_samp.iloc[row, 4]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum RT (min)'] = total_samp.iloc[
                                                          row, 1]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum Precursor m/z'] = total_samp.iloc[
                                                               row, 0]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum sn'] = total_samp.iloc[row, 2]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum score'] = total_samp.iloc[row, 3]
            else:
                # if index changes, moves to next sample_name column
                col_index += 1
                overlap = np.where(((alignment_df.iloc[:, 1] - RT_error)
                                   <= total_samp.iloc[row, 1]) &
                                   (total_samp.iloc[row, 1] <=
                                   (alignment_df.iloc[:, 1] + RT_error))
                                   & ((alignment_df.iloc[:, 0] - mz_error)
                                      <= total_samp.iloc[row, 0])
                                   & (total_samp.iloc[row, 0] <=
                                       (alignment_df.iloc[:, 0] + mz_error)))
                if len(overlap[0]) > 0:
                    overlap = [overlap][0][0][0]
                    if (alignment_df.iloc[overlap, col_index]) > 0:
                        alignment_df.at[overlap, alignment_df.columns[
                                         col_index]] += total_samp.iloc[row, 4]
                        # Averages the mean if there is a repeat intensity
                        alignment_df.at[overlap, alignment_df.columns[
                                         col_index]] /= 2
                    else:
                        alignment_df.at[overlap, alignment_df.columns[
                                             col_index]] = total_samp.iloc[
                                                           row, 4]
                        # adds rt, m/z, score, sn to aligned point
                        alignment_df.at[overlap,
                                        'Sum RT (min)'] += (total_samp.
                                                            iloc[row, 1])
                        alignment_df.at[overlap,
                                        'Sum Precursor m/z'] += (total_samp.
                                                                 iloc[row, 0])
                        alignment_df.at[overlap, 'Sum sn'] += (total_samp.
                                                               iloc[row, 2])
                        alignment_df.at[overlap, 'Sum score'] += (total_samp.
                                                                  iloc[row, 3])
                else:
                    # appends as a new reference if nothing is found
                    alignment_df.loc[row_max] = total_samp.iloc[row, 0:4]
                    alignment_df.at[alignment_df.index[
                                   -1], alignment_df.columns[
                                        col_index]] = total_samp.iloc[row, 4]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum RT (min)'] = total_samp.iloc[
                                                          row, 1]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum Precursor m/z'] = total_samp.iloc[
                                                               row, 0]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum sn'] = total_samp.iloc[row, 2]
                    alignment_df.at[alignment_df.index[
                                   -1], 'Sum score'] = total_samp.iloc[row, 3]
        else:  # for last row of dataframe
            overlap = np.where(((alignment_df.iloc[:, 1] - RT_error)
                               <= total_samp.iloc[row, 1]) &
                               (total_samp.iloc[row, 1] <=
                                (alignment_df.iloc[:, 1] + RT_error))
                               & ((alignment_df.iloc[:, 0] - mz_error)
                                   <= total_samp.iloc[row, 0])
                               & (total_samp.iloc[row, 0]
                                   <= (alignment_df.iloc[:, 0] + mz_error)))
            if len(overlap[0]) > 0:
                overlap = [overlap][0][0][0]
                if (alignment_df.iloc[overlap, col_index]) > 0:
                    alignment_df.at[overlap, alignment_df.columns[
                                         col_index]] += total_samp.iloc[row, 4]
                    # Averages the mean if there is a repeat intensity
                    alignment_df.at[overlap, alignment_df.columns[
                                         col_index]] /= 2
                else:
                    # sets intensity value where realignment was identified
                    alignment_df.at[overlap, alignment_df.columns[
                                         col_index]] = total_samp.iloc[row, 4]
                    # adds rt, m/z, score, sn to aligned point
                    alignment_df.at[overlap, 'Sum RT (min)'] += (total_samp.
                                                                 iloc[row, 1])
                    alignment_df.at[overlap,
                                    'Sum Precursor m/z'] += (total_samp.
                                                             iloc[row, 0])
                    alignment_df.at[overlap, 'Sum sn'] += total_samp.iloc[
                                                           row, 2]
                    alignment_df.at[overlap, 'Sum score'] += total_samp.iloc[
                                                              row, 3]
            else:
                # appends as a new reference if nothing is found
                alignment_df.loc[row_max] = total_samp.iloc[row, 0:4]
                alignment_df.at[alignment_df.index[
                               -1], alignment_df.columns[
                                    col_index]] = total_samp.iloc[row, 4]
                alignment_df.at[alignment_df.index[
                               -1], 'Sum RT (min)'] = total_samp.iloc[row, 1]
                alignment_df.at[alignment_df.index[
                               -1], 'Sum Precursor m/z'] = total_samp.iloc[
                                                           row, 0]
                alignment_df.at[alignment_df.index[
                               -1], 'Sum sn'] = total_samp.iloc[row, 2]
                alignment_df.at[alignment_df.index[
                               -1], 'Sum score'] = total_samp.iloc[row, 3]
    alignment_df.rename(columns={'rt': 'Average RT (min)'}, inplace=True)
    alignment_df.rename(columns={'m/z': 'Average m/z'}, inplace=True)
    alignment_df.rename(columns={'sn': 'Average sn'}, inplace=True)
    alignment_df.rename(columns={'score': 'Average score'}, inplace=True)
    # Replace all NaN elements with 0
    alignment_df = alignment_df.fillna(0)
    # final dataframe sorted by m/z
    alignment_df = alignment_df.sort_values(by='Average m/z',
                                            ignore_index=True)
    # empty list to hold invalid peaks
    invalid = []
    for rows in range(len(alignment_df)):
        # Calculating the averages after the iterations
        # requires count of nonzero count to calculate the mean properly
        count = np.count_nonzero(alignment_df.iloc[rows, 8:])
        if count > 0:
            alignment_df.at[rows, 'Average RT (min)'] = (alignment_df.loc[rows,
                                                         'Sum RT (min)']/count)
            alignment_df.at[rows, 'Average m/z'] = (alignment_df.loc[rows,
                                                    'Sum Precursor m/z']/count)
            alignment_df.at[rows, 'Average sn'] = (alignment_df.loc[rows,
                                                   'Sum sn']/count)
            alignment_df.at[rows, 'Average score'] = (alignment_df.loc[rows,
                                                      'Sum score']/count)
        else:
            invalid.append(rows)
            # print("Identified invalid peak(s)")  # For detection error
    if len(invalid) > 0:
        alignment_df = alignment_df.drop(invalid)
        alignment_df = alignment_df.reset_index(drop=True)
        print("Removed invalid peak(s) from reference")
    else:
        pass
    # Drop columns to collect sums for averaging
    alignment_df = alignment_df.drop(columns=[
                   'Sum RT (min)', 'Sum Precursor m/z',
                   'Sum sn', 'Sum score'])
    # rounds numbers in the first two columns
    alignment_df = alignment_df.round({'Average RT (min)': 3,
                                       'Average m/z': 5})
    print("Alignment done!")
    # converts file for saving
    if export_name == None:
        print("Result didn't exported!")
        pass
    else:
        alignment_df.to_csv(export_name,
                            header=True, index=False)
    print("Result saved to " + str(export_name))
    return alignment_df


def mss_process(path, export_name, align_mz_err=0.015, align_rt_err=0.1,
              remove_noise=True, thres_noise=1000,
              err_ppm=10, enable_score=True, mz_c_thres=5, peak_base=0.005,
              peakutils_thres=0.02, min_d=1, rt_window=1.5,
              peak_area_thres=1e4, min_scan=5, max_scan=50,
              max_peak=5):
    print('Reading files...')
    batch_scan, name_list = batch_scans(path, remove_noise=remove_noise,
                                            thres_noise=thres_noise)
    print('Processing peak list...')
    d_peak = []
    for i in range(len(batch_scan)):
        print('Processing', str(int(i + 1)), 'out of ', len(batch_scan), 'file')
        d_result = peak_list(batch_scan[i], err_ppm=err_ppm,
                            enable_score=enable_score, mz_c_thres=mz_c_thres,
                            peak_base=peak_base, peakutils_thres=peakutils_thres,
                            min_d=min_d, rt_window=rt_window,
                            peak_area_thres=peak_area_thres, min_scan=min_scan,
                            max_scan=max_scan, max_peak=max_peak)
        d_peak.append(d_result)

    d_align = realignment(d_peak, export_name, name_list, align_rt_err=align_rt_err,
                        align_mz_err=align_mz_err)
    
    print('Finished!')
    return d_align