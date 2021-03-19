import sys
# import mss
sys.path.append('../')
from mss import mssmain as msm
from mss import align
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle


def s_stack(d_batch):
    """This function compiles all samples files into one dataframe
    for analysis. It takes in files that are of the .txt type."""
    all_samples = []
    # Calculating the number of files in path
    num_files = len(d_batch)
    for i in range(num_files):
        sample_name = d_batch[i]
        sample_df = pd.DataFrame(sample_name, columns=['m/z',
                                                       'rt', 'sn',
                                                       'score', 'peak area'],
                                 dtype=np.float32)
        sample_df = sample_df.sample(frac=1).reset_index(drop=True)
        all_samples += [sample_df]
    # Combining all the dataframes into one
    tot_samples = pd.concat(all_samples)
    all_samples.clear()
    # Cleaning up the data before processing
    tot_samples.loc[tot_samples.sn == 0, 'sn'] = 5000.0  # 5000 as max
    tot_samples.loc[tot_samples.score == 3.0, 'score'] = 0.1
    tot_samples.loc[tot_samples.score == 2.0, 'score'] = 0.6
    print("Process completed!")
    return tot_samples


def mss_align(d_batch, export_name, name_list, RT_error, mz_error):
    np.warnings.filterwarnings('ignore')  # ignore warning from sn inf value
    tot_samp = s_stack(d_batch)
    alignment_df = tot_samp.iloc[:1, :4].copy()
    alignment_df[['s_rt', 's_mz', 's_sn', 's_scr']] = [0.0, 0.0, 0.0, 0.0]
    alignment_df = alignment_df.astype({'s_rt': 'float32',
                                        's_mz': 'float32', 's_sn': 'float32',
                                        's_scr': 'float32'})
    # col_index to track where to add sample_name columns
    col_index = 8
    name_list = [i[:-5] for i in name_list]
    alignment_df[name_list] = 0.0
    alignment_df[name_list] = alignment_df[name_list].astype('float32')
    print("Alignment beginning..")

    def value_set(df, loc, row, col_index):
        # sets intensity value where realignment was identified
        c = df.columns[col_index]
        df.at[loc, c] = tot_samp.iloc[row, 4]
        df.iloc[loc] = np.nan_to_num(df.iloc[loc])
        df.at[loc, 's_mz'] += tot_samp.iloc[row, 0]
        df.at[loc, 's_rt'] += tot_samp.iloc[row, 1]
        df.at[loc, 's_sn'] += tot_samp.iloc[row, 2]
        df.at[loc, 's_scr'] += tot_samp.iloc[row, 3]
        return

    for row in tqdm(range(len(tot_samp))):
        if alignment_df.isnull().values.any() is True:
            alignment_df = alignment_df.fillna(0)
        else:
            pass
        row_max = len(alignment_df)
        if row < range(len(tot_samp))[-1]:
            # Loops until a change in sample_name occurs by index
            if tot_samp.index[row] < tot_samp.index[row + 1]:
                # checks to see if row exists in reference or not
                overlap = np.where((abs(alignment_df.iloc[:, 1] -
                                        tot_samp.iloc[row, 1])
                                    <= RT_error)
                                   & (abs(alignment_df.iloc[:, 0] -
                                          tot_samp.iloc[row, 0])
                                      <= mz_error))
                if len(overlap[0]) > 0:
                    overlap = [overlap][0][0][0]
                    # checks for any duplicates
                    if (alignment_df.iloc[overlap, col_index]) > 0:
                        alignment_df.at[overlap, alignment_df.columns[
                                        col_index]] += tot_samp.iloc[row, 4]
                        # Averages the mean if there is a repeat intensity
                        alignment_df.at[overlap, alignment_df.columns[
                            col_index]] /= 2
                    else:
                        # sets intensity value where realignment was identified
                        value_set(alignment_df, overlap, row, col_index)
                else:
                    # appends as a new reference if nothing is found
                    alignment_df.loc[row_max] = tot_samp.iloc[row, 0:4]
                    value_set(alignment_df,
                              alignment_df.index[-1], row, col_index)
            else:
                # if index changes, moves to next sample_name column
                col_index += 1
                overlap = np.where((abs(alignment_df.iloc[:, 1] -
                                        tot_samp.iloc[row, 1])
                                    <= RT_error)
                                   & (abs(alignment_df.iloc[:, 0] -
                                          tot_samp.iloc[row, 0])
                                      <= mz_error))
                if len(overlap[0]) > 0:
                    overlap = [overlap][0][0][0]
                    if (alignment_df.iloc[overlap, col_index]) > 0:
                        alignment_df.at[overlap, alignment_df.columns[
                            col_index]] += tot_samp.iloc[row, 4]
                        # Averages the mean if there is a repeat intensity
                        alignment_df.at[overlap, alignment_df.columns[
                            col_index]] /= 2
                    else:
                        value_set(alignment_df, overlap, row, col_index)
                else:
                    # appends as a new reference if nothing is found
                    alignment_df.loc[row_max] = tot_samp.iloc[row, 0:4]
                    value_set(alignment_df,
                              alignment_df.index[-1], row, col_index)
        else:  # for last row of dataframe
            overlap = np.where(((alignment_df.iloc[:, 1] - RT_error)
                                <= tot_samp.iloc[row, 1]) &
                               (tot_samp.iloc[row, 1] <=
                                (alignment_df.iloc[:, 1] + RT_error))
                               & ((alignment_df.iloc[:, 0] - mz_error)
                                   <= tot_samp.iloc[row, 0])
                               & (tot_samp.iloc[row, 0]
                                   <= (alignment_df.iloc[:, 0] + mz_error)))
            if len(overlap[0]) > 0:
                overlap = [overlap][0][0][0]
                if (alignment_df.iloc[overlap, col_index]) > 0:
                    alignment_df.at[overlap, alignment_df.columns[
                        col_index]] += tot_samp.iloc[row, 4]
                    # Averages the mean if there is a repeat intensity
                    alignment_df.at[overlap, alignment_df.columns[
                        col_index]] /= 2
                else:
                    value_set(alignment_df, overlap, row, col_index)
            else:
                # appends as a new reference if nothing is found
                value_set(alignment_df, alignment_df.index[-1], row, col_index)
    alignment_df.rename(columns={'rt': 'Average_rt',
                                 'm/z': 'Average m/z',
                                 'sn': 'Average sn',
                                 'score': 'Average score'}, inplace=True)
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
            alignment_df.at[rows, 'Average_rt'] = (
                alignment_df.loc[rows, 's_rt']/count)
            alignment_df.at[rows, 'Average m/z'] = (
                alignment_df.loc[rows, 's_mz']/count)
            alignment_df.at[rows, 'Average sn'] = (
                alignment_df.loc[rows, 's_sn']/count)
            alignment_df.at[rows, 'Average score'] = (
                alignment_df.loc[rows, 's_scr']/count)
        else:
            invalid.append(rows)
            print("Identified invalid peak(s)")  # For detection error
    if len(invalid) > 0:
        alignment_df = alignment_df.drop(invalid)
        alignment_df = alignment_df.reset_index(drop=True)
    else:
        pass
    # Drop columns to collect sums for averaging
    alignment_df = alignment_df.drop(columns=[
        's_rt', 's_mz',
        's_sn', 's_scr'])
    # rounds numbers in the first two columns
    alignment_df = alignment_df.round({'Average_rt': 3,
                                       'Average m/z': 5})
    print("Alignment done!")
    # converts file for saving
    if export_name is None:
        print("Result didn't exported!")
        pass
    else:
        alignment_df.to_csv(export_name,
                            header=True, index=False)
    print("Result saved to " + str(export_name))
    return alignment_df


print('Reading files...')
path = 'D:/UW/directproject/example_data/'
batch_scan, name_list = msm.batch_scans(path, remove_noise=True, thres_noise=20000)
print('Processing peak list...')
d_peak = []
for i in range(len(batch_scan)):
    print('Processing', str(int(i + 1)),
            'out of ', len(batch_scan), 'file')
    d_result = msm.peak_list(batch_scan[i], enable_score=False)
    d_peak.append(d_result)

shape=[]
for i in range(100):
    d_align = mss_align(d_peak, export_name=None, name_list=name_list, RT_error=0.1, mz_error=0.015)
    shape.append(d_align.shape[0])
    print(i)

with open('1000loop.shape', 'wb') as fp:
	pickle.dump(shape, fp)

print('Finished!')
