import numpy as np
import h5py as h5
from scipy.signal import medfilt


def mask_freq_RFI(vis, sig=2., size=21):

    vis_tmean = np.ma.mean(vis, axis=0)
    vis_fsmooth = vis_tmean - medfilt(vis_tmean, [size, 1, 1])
    f_var  = np.ma.var(vis_fsmooth, axis=0)
    f_mean = np.ma.mean(vis_fsmooth, axis=0)
    f_mask  = vis_fsmooth >  (f_mean + sig * f_var)[None, ...]
    f_mask += vis_fsmooth <  (f_mean - sig * f_var)[None, ...]

    return f_mask


def load_data_list(data_path, data_file_list, freq_sel = [None, None], bad_freq=[]):

    data = None
    for data_file in data_file_list:
        data = load_data(data_path, data_file, data=data,
                         freq_sel = freq_sel, bad_freq=bad_freq)

    return data

def load_data(data_base, data_file, data=None, freq_sel = [None, None], bad_freq=[]):

    with h5.File(data_base + data_file, 'r') as fh: 
        freqstart = fh.attrs['freqstart']
        freqstep  = fh.attrs['freqstep']
        freqn     = fh.attrs['nfreq']
        freq = np.arange(freqn) * freqstep + freqstart
        
        ants = fh['blorder'][:]
        
        freq = freq[slice(*freq_sel)]
    
        vis = fh['vis'][:, slice(*freq_sel), ...].real
    
        timestart = fh.attrs['sec1970']
        timestep  = fh.attrs['inttime']
        timen     = fh['vis'].shape[0]
        time = np.arange(timen) * timestep + timestart
        
        ra = fh['ra'][:]
        dec= fh['dec'][:]
    
    vis = np.ma.array(vis)
    vis.mask = np.zeros(vis.shape).astype('bool')
    for _bad_freq in bad_freq:
        vis.mask[:, slice(*_bad_freq), ...] = True
    f_mask = mask_freq_RFI(vis)
    vis.mask += f_mask[None, ...]

    if data is None:
        data = {}
        data['vis']  = vis
        data['freq'] = freq
        data['time'] = time
        data['ants'] = ants
        data['ra'] = ra
        data['dec'] = dec
    else:
        data['vis']  = np.ma.concatenate([data['vis'], vis],   axis=0)
        data['time'] = np.ma.concatenate([data['time'], time], axis=0)
        data['ra']   = np.ma.concatenate([data['ra'], ra],     axis=0)
        data['dec']  = np.ma.concatenate([data['dec'], ra],    axis=0)

    return data


