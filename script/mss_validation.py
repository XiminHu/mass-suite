import sys
sys.path.append('../mss')
import matplotlib.pyplot as plt
import pandas as pd
import visreader as mvis
import mssmain as mss
import numpy as np

sample = input('sample number:')
path = 'D:/UW/mssvalidation/20201119_ENTACT_validation_data/20201117_ENTACT_' + str(sample) + '_1.mzML'
scans = mss.get_scans(path, ms_all=False, ms_lv=1)
#noise removal
mss.noise_removal(scans, 2000)

df_path = 'D:/UW/mssvalidation/20201119_ENTACT_validation_data/' + str(sample) + '_minscan5.csv'
df = pd.read_csv(df_path)
df.drop(columns=['Unnamed: 0'], inplace=True)
d_val = pd.read_excel('D:/UW/mssvalidation/20201119_ENTACT_validation_data/ENTACT_information/NTA_DataReporting_UW_CUW_unblinded_internal.xlsx', sheet_name=sample)

d_val = d_val[d_val['Ion Polarity'] == 'Positive']

RT_error = 2
mz_error = 0.015
for i,row in d_val.iterrows():
    print(i, 'out of ', d_val.index[-1])
    mz_ref = row['Reference mass'] + 1.0079
    rt_ref = row['RT']
    print('ref mz: ', mz_ref, 'ref RT: ', rt_ref)
    overlap = np.where(((abs(df['m/z'] - mz_ref) <= mz_error)) & ((abs(df['rt'] - rt_ref)) <= RT_error))
    print(df.loc[overlap])
    mvis.integration_plot(scans, mz_ref, 10,fig_w=10,fig_h=3)
    plt.show()
    d_val.at[i, 'match'] = input('result, 1 for match, 2 for dont know, 3 for no match:\r')
    print(' ')
    print(' ')

d_val.to_csv('D:/UW/mssvalidation/20201119_ENTACT_validation_data/export/ENTACT_match/' + str(sample) + '.csv')