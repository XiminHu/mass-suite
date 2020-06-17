import sys
sys.path.append('../mss')
import matplotlib.pyplot as plt

import visreader as mvis
import mssmain as mss
import pandas as pd
import numpy as np
from tqdm import tqdm

def mz_selector(scans, mz_list, export_name):
    sample_list = []
    for mz in tqdm(mz_list):
        rt, i = ms_chromatogram_list(scans, mz, 20)
        count = 0
        for ints in i:
            if ints >= 5000:
                count += 1
                if count == 7:
                    sample_list.append([mz, i])
                    break
            else:
                count = 0
                continue
                
    d_sample = pd.DataFrame(sample_list)
    d_rt = pd.DataFrame(rt)
    
    writer = pd.ExcelWriter(export_name, engine='xlsxwriter')
    d_sample.to_excel(writer, sheet_name='intensity')
    d_rt.to_excel(writer, sheet_name='rt')
    writer.save()
    
    return


path = input('please input the path of mzml files:')
batch_scan, file_list = mss.batch_scans(path)
for scans in batch_scan:
    mss.noise_removal(scans, 5000)

output = input('output path:')
for i in range(len(batch_scan)):
    print(str(i+1), 'out of', len(batch_scan))
    mz_list = mz_gen(batch_scan[i])
    output_name = output + str(i + 1) + '.xlsx'
    mz_selector(batch_scan[i], mz_list,output_name)
print('finished!')