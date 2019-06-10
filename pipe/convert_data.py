from FAST_utils.src import data_conv as dc

data_path   = '/scratch/users/ycli/.test/'
output_path = data_path + 'test.h5'
data_file   = 'SDSS_N_2a_arcdrfit-M%02d_W_%04d.fits'
dec0 = 25.65294
beam_list   = [1, 2]
block_list  = [1, 2]

dc.convert_to_tl(data_path, data_file, output_path,
        dec0 = 25.65294, feed_rotation=0., 
        beam_list=beam_list, block_list=block_list, 
        degrade_freq_resol=16)
