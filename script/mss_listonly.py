import sys
sys.path.append('../mss')
import matplotlib.pyplot as plt
import visreader as mvis
import mssmain as mss

path = input('Please input the mzml file path:')
noise_thres = int(input('Please input the noise threshold for ms1 spectrum:'))

print('Reading mzml files...')
batch_scan, file_list = mss.batch_scans(path, True, noise_thres)

output_dir = input('Please input the export path')
print('Processing peak list...')
for i in range(len(batch_scan)):
	print('Processing', str(int(i+1)), 'out of ', len(batch_scan), 'file')
	d_result = mss.peak_list(batch_scan[i], 20)
	output_name = file_list[i][:-5] + '.csv'
	output_path = output_dir + output_name
	d_result.to_csv(output_path)

print('Finished!')