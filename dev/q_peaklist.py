import sys
# import mss
sys.path.append('../')
from mss import mssmain as msm
from mss import align

path = input('filename to process:')
scans = msm.get_scans(path, ms_all=False, ms_lv=1)
#noise removal
msm.noise_removal(scans, 2000)
d_op = msm.peak_list(scans, 10, enable_score=True,peak_base=0.001,peak_area_thres=0)
output_path = input('path and filename to export:')
d_op.to_csv(output_path)
