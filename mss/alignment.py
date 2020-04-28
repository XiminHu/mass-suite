import pandas as pd
import numpy as np
from tqdm import tqdm
import os
'''
coordination: as most of the position argument is now using the relative position rather than label, 
later coordination need to make sure new values won't interrupt position of existing values
'''

def stack (batch_path):
    """This function compiles all samples files into one dataframe for analysis."""
    all_samples = []
    num_files = len([f for f in os.listdir(batch_path)if os.path.isfile(os.path.join(batch_path, f))])
    for files in tqdm(range(num_files)):
        sample = os.listdir(batch_path)[files]
        sample_df = pd.read_csv(batch_path + sample, sep="\\t", usecols=['RT (min)', 'Precursor m/z', 'Area'], dtype=np.float32, engine='python')
        all_samples += [sample_df]
    total_samples = pd.concat(all_samples)
    all_samples.clear()
    return total_samples

def realignment(batch_path, batch_name):
#def realignment(batch_path, batch_name, rt_error=0.05, mz_error=0.015)
#maybe add some print option to make the script more explicit

    """This function works by using one .txt file as a reference in which
    other files realigned to in terms of precursor and RT. """
    RT_error = 0.05 # units of minutes, can be adjusted
    alignment_df = pd.DataFrame()
    standard_sample = os.listdir(batch_path)[0] # first sample
    # reads .txt file into a dataframe
    standard_df = pd.read_csv(batch_path + standard_sample, sep="\\t", usecols=['RT (min)', 'Precursor m/z'], dtype= np.float32, engine='python')
    alignment_df['RT (min)'] = standard_df.iloc[:,0]
    alignment_df['Precursor m/z'] = standard_df.iloc[:,1]
    alignment_df['Sum RT (min)'] = 0.0
    alignment_df['Sum Precursor m/z'] = 0.0
    alignment_df['Sum RT (min)'] = alignment_df['Sum RT (min)'].astype('float32')
    alignment_df['Sum Precursor m/z'] = alignment_df['Sum Precursor m/z'].astype('float32')
    mz_error = 0.015
#col_index will be coordinated later to bring in more values if needed
    col_index = 4
#maybe take out the progress bar in stack function and instead use
#print('reading in files')
    total_samp = stack(batch_path)
<<<<<<< HEAD
#maybe leave a backddor in stack function to get the num_files
#in stack, instead of return the df, try return num_files, df
#Thus you will have the function as total_samp, num_files = stack(batch_path)
    num_files = len([f for f in os.listdir(batch_path)if os.path.isfile(os.path.join(batch_path, f))]) 
=======
    num_files = len([f for f in os.listdir(batch_path)if os.path.isfile(os.path.join(batch_path, f))])
>>>>>>> 22819954d4ffc78f45a5d9a4d10e42c5855f4016
    for files in range(num_files):
#sample string trim is needed to get rid of the .txt/.csv extension
        sample = os.listdir(batch_path)[files] #Create columns using the file name
        alignment_df[sample] = 0.0
        alignment_df[sample] = alignment_df[sample].astype('float32')
    for row in tqdm(range(len(total_samp))):
        row_max = len(alignment_df)
        if row < range(len(total_samp))[-1]:
            if total_samp.index[row] < total_samp.index[row + 1]:
                overlap = np.where(((alignment_df.iloc[:, 0] - RT_error) <= total_samp.iloc[row, 0]) &\
                                   (total_samp.iloc[row, 0] <= (alignment_df.iloc[:, 0] + RT_error))\
                                     & ((alignment_df.iloc[:, 1] - mz_error) <= total_samp.iloc[row, 1])\
                                     & (total_samp.iloc[row, 1] <= (alignment_df.iloc[:, 1] + mz_error)))
                if len(overlap[0])>0:
                    overlap = [overlap][0][0][0]
                    # sets intensity value into slot where realignment was identified
                    alignment_df.at[overlap, alignment_df.columns[col_index]] = total_samp.iloc[row, 2] #Merged the latest one, potentially the closest one?
                    # adds RT time to aligned point
                    alignment_df.at[overlap, 'Sum RT (min)'] += total_samp.iloc[row, 0]
                    # adds m/z value to aligned point
                    alignment_df.at[overlap, 'Sum Precursor m/z'] += total_samp.iloc[row, 1]
                else:
                    # appends as a new reference if nothing is found
                    alignment_df.loc[row_max] = total_samp.iloc[row, 0:2]
                    alignment_df.at[alignment_df.index[-1], alignment_df.columns[col_index]] = total_samp.iloc[row, 2]
                    alignment_df.at[alignment_df.index[-1], 'Sum RT (min)'] = total_samp.iloc[row, 0]
                    alignment_df.at[alignment_df.index[-1], 'Sum Precursor m/z'] = total_samp.iloc[row, 1]
            else:
                col_index += 1
                overlap = np.where(((alignment_df.iloc[:, 0] - RT_error) <= total_samp.iloc[row, 0]) &\
                                   (total_samp.iloc[row, 0] <= (alignment_df.iloc[:, 0] + RT_error))\
                                     & ((alignment_df.iloc[:, 1] - mz_error) <= total_samp.iloc[row, 1])\
                                     & (total_samp.iloc[row, 1] <= (alignment_df.iloc[:, 1] + mz_error)))
                if len(overlap[0])>0:
                    overlap = [overlap][0][0][0]
                    # sets intensity value into slot where realignment was identified
                    alignment_df.at[overlap, alignment_df.columns[col_index]] = total_samp.iloc[row, 2]
                    # adds RT time to aligned point
                    alignment_df.at[overlap, 'Sum RT (min)'] += total_samp.iloc[row, 0]
                    # adds m/z value to aligned point
                    alignment_df.at[overlap, 'Sum Precursor m/z'] += total_samp.iloc[row, 1]
                else:
                    # appends as a new reference if nothing is found
                    alignment_df.loc[row_max] = total_samp.iloc[row, 0:2]
                    alignment_df.at[alignment_df.index[-1], alignment_df.columns[col_index]] = total_samp.iloc[row, 2]
                    alignment_df.at[alignment_df.index[-1], 'Sum RT (min)'] = total_samp.iloc[row, 0]
                    alignment_df.at[alignment_df.index[-1], 'Sum Precursor m/z'] = total_samp.iloc[row, 1]
        else:
            overlap = np.where(((alignment_df.iloc[:, 0] - RT_error) <= total_samp.iloc[row, 0]) &\
                                (total_samp.iloc[row, 0] <= (alignment_df.iloc[:, 0] + RT_error))\
                                 & ((alignment_df.iloc[:, 1] - mz_error) <= total_samp.iloc[row, 1])\
                                 & (total_samp.iloc[row, 1] <= (alignment_df.iloc[:, 1] + mz_error)))
            if len(overlap[0])>0:
                overlap = [overlap][0][0][0]
                # sets intensity value into slot where realignment was identified
                alignment_df.at[overlap, alignment_df.columns[col_index]] = total_samp.iloc[row, 2]
                # adds RT time to aligned point
                alignment_df.at[overlap, 'Sum RT (min)'] += total_samp.iloc[row, 0]
                # adds m/z value to aligned point
                alignment_df.at[overlap, 'Sum Precursor m/z'] += total_samp.iloc[row, 1]
            else:
                # appends as a new reference if nothing is found
                alignment_df.loc[row_max] = total_samp.iloc[row, 0:2]
                alignment_df.at[alignment_df.index[-1], alignment_df.columns[col_index]] = total_samp.iloc[row, 2]
                alignment_df.at[alignment_df.index[-1], 'Sum RT (min)'] = total_samp.iloc[row, 0]
                alignment_df.at[alignment_df.index[-1], 'Sum Precursor m/z'] = total_samp.iloc[row, 1]
    alignment_df.rename(columns = {'RT (min)':'Average RT (min)'}, inplace = True)
    alignment_df.rename(columns = {'Precursor m/z':'Average m/z'}, inplace = True)
    # Replace all NaN elements with 0
    alignment_df = alignment_df.fillna(0)
    for rows in range(len(alignment_df)):
        # Calculating the averages after the iterations
        # requires count of nonzero count to calculate the mean properly

#Also note here for later update since more value could be inserted, so instead of .iloc[,4:] it could be different
        count = np.count_nonzero(alignment_df.iloc[rows, 4:])
        if count > 0:
            alignment_df.at[rows, 'Average RT (min)'] = alignment_df.loc[rows, 'Sum RT (min)']/count
            alignment_df.at[rows, 'Average m/z'] = alignment_df.loc[rows, 'Sum Precursor m/z']/count

#Try throw a exception/text here if count==0, since it could be a debug need or detection on failure of oritinal files
        else:
            pass
    # Drop columns to collect sums for averaging
    alignment_df = alignment_df.drop(columns=['Sum RT (min)', 'Sum Precursor m/z'])
    # Final sort by m/z
#got error: sort_values() got an unexpected keyword argument 'ignore_index', pandas ver 0.25.1
#after updated to pandas 1.0.3 works
    alignment_df = alignment_df.sort_values(by='Average m/z', ignore_index=True)
    # converts file for saving
#Leave out a backdoor for user to define the output path
#alignment_df.to_csv(path + batch_name + '.csv', header=True, index=False)
    alignment_df.to_csv(batch_name + '.csv', header=True, index=False)
    return alignment_df
