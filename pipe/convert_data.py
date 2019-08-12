import numpy as np
from FAST_utils.src import data_conv as dc
import os

data_sect   = '/SDSS_N_2a/20190527/'
#data_sect   = '/SDSS_N_2a/20190529/'
data_file_tmp = 'SDSS_N_2a_arcdrfit'

data_path   = '/data31/2019a-053-SC/'
data_file   = data_sect + data_file_tmp + '-M%02d_W_%04d.fits'

block_st = 11
block_ed = 20

output_path = '/home/wangyg/data/' + data_sect
if not os.path.exists(output_path):
    os.makedirs(output_path)
output_file = data_file_tmp + '%04d-%04d.h5'%(block_st, block_ed)

#dec0 = 26.15294 
# for data on 27
alt=89.93462193661559
az=180.
beam_list   = np.arange(19) + 1
block_list  = np.arange(block_st, block_ed + 1)

dc.convert_to_tl(data_path, data_file, 
        output_path + output_file,
        alt=alt, az=az, feed_rotation=0., 
        beam_list=beam_list, block_list=block_list, 
        degrade_freq_resol=128)
        #degrade_freq_resol=16)
