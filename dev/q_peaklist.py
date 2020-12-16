import sys
#import mss
print('__file__={0:<35} | __name__={1:<20} | __package__={2:<20}'.format(__file__,__name__,str(__package__)))
sys.path.append('../mss')
import visreader as mvis
import mssmain as mss
import pandas as pd

path = input('filename to process:')
scans = mss.get_scans(path, ms_all=False, ms_lv=1)
#noise removal
mss.noise_removal(scans, 2000)
d_op = mss.peak_list(scans, 10, enable_score=True,peak_base=0.001,peak_area_thres=0)
output_path = input('path and filename to export:')
d_op.to_csv(output_path)
