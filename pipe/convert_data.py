import numpy as np
from FAST_utils.src import data_conv as dc

data_path   = '/data31/2019a-053-SC/SDSS_N_2a/20190527/'
data_file_tmp = 'SDSS_N_2a_arcdrfit'
output_path = '/home/wangyg/tmp/' + data_file_tmp + '%04d-%04d.h5'
data_file   = data_file_tmp + '-M%02d_W_%04d.fits'
#dec0 = 25.65294
dec0 = 26.15294 
beam_list   = np.arange(19) + 1
block_list  = np.arange(10) + 1 + 10

dc.convert_to_tl(data_path, data_file, 
        output_path%(block_list[0], block_list[-1]),
        dec0 = 25.65294, feed_rotation=0., 
        beam_list=beam_list, block_list=block_list, 
        degrade_freq_resol=1024)
