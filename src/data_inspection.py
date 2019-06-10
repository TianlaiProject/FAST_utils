#import pyfits
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator,ScalarFormatter


import h5py

from astropy import units as u
from astropy.time import Time

from datetime import datetime, timedelta
import matplotlib.dates as mdates
from scipy.signal import medfilt
from scipy.ndimage import median_filter
from scipy.signal import lombscargle
import copy

_c_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "0.2", "#d62728", "#9467bd",
                  "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", 'w']
_l_list = ['s-', 'o--', '.--', 's--', 'o-']


#class DATA_Inspection(object):


class DATA_Inspection_old(object):
    
    def get_polymean(self, deg=6):
        
        freq = cls.freq
        #data = self.data
        #data = data.reshape((-1, ) + data.shape[-2:])
        spec_mean = np.ma.mean(data, axis=0)
        spec_poly = np.zeros_like(spec_mean)
        poly_fit  = np.polyfit(freq, spec_mean.T, deg, rcond=1.e-10)
        for ii in range(deg + 1):
            spec_poly += freq[None, :] ** (deg - ii) * poly_fit[ii][:, None]
        #poly_fit  = 10.**poly_fit
        
        return spec_poly
        #return data - spec_poly[None, :, :], spec_poly
        
    def get_bandpass(self, data, kernel=101):
        
        spec_mean = np.ma.mean(data, axis=0)
        time_mask = np.all(data.mask, axis=(1, 2))
        freq_mask = np.sum(data.mask[~time_mask, ...], axis=(0,1)).astype('bool')
        #print freq_mask.shape
        #print freq_mask
        
        #print spec_mean.shape
        spec_mean[:, ~freq_mask] = median_filter(spec_mean[:, ~freq_mask], size = (1, kernel))
        spec_mean = np.ma.array(spec_mean)
        spec_mean.mask[:, freq_mask] = True
        
        return spec_mean

    
    def rfi_flagging(self, sigma=6., bad_time_list=[], bad_freq_list=[],
                     loop = 2, mask_width=10):
        
        data = np.ma.array(self.data)
        mask = np.zeros(data.shape).astype('bool') + self.mask
        
        data = data.reshape((-1, ) + data.shape[-2:])
        mask = mask.reshape(data.shape)
        #bdps = self.get_bandpass(data, kernel=301, freq_mask=None)
        #data = data / bdps[None, :, :]
        
        
        for bad_time in bad_time_list:
            print "mask time slice of ", bad_time
            mask[slice(*bad_time)] = True
        for bad_freq in bad_freq_list:
            print "mask freq slice of ", bad_freq
            mask[:,:,slice(*bad_freq)] = True
        print np.any(mask)
        
        data.mask = mask
        
        freq_mask = np.zeros(self.freq.shape).astype('bool')
        #print freq_mask.shape
        mask_rate_raw = np.sum(freq_mask)
        
        bdps = self.get_bandpass(data, kernel=301)
        spec_mean = np.ma.mean(data / bdps[None, :, :], axis=0)
        for ll in range(loop):

            spec_mean_mean = np.ma.mean(spec_mean[:,~freq_mask], axis=-1)
            spec_mean_std  = np.ma.std(spec_mean[:,~freq_mask],  axis=-1)
            print spec_mean_mean, spec_mean_std
        
            _freq_mask  = spec_mean[0, :] > (spec_mean_mean[0] + sigma * spec_mean_std[0])
            _freq_mask += spec_mean[0, :] < (spec_mean_mean[0] - sigma * spec_mean_std[0])
            _freq_mask += spec_mean[1, :] > (spec_mean_mean[1] + sigma * spec_mean_std[1])
            _freq_mask += spec_mean[1, :] < (spec_mean_mean[1] - sigma * spec_mean_std[1])
            
            _freq_mask += np.roll(_freq_mask, mask_width/2) + np.roll(_freq_mask, -mask_width/2)
            
            freq_mask += _freq_mask
            
            mask_rate = np.sum(freq_mask)
            print "Mask Add %d"%(mask_rate - mask_rate_raw)
            mask_rate_raw = mask_rate
        
        mask += freq_mask[None, None, :]
        data.mask = mask
        
        #bdps = self.get_bandpass(data, kernel=51)
        #self.plot_spec(data = data/bdps[None, :, :], nblock=None)
        #self.plot_spec(data = data, nblock=None)
        #_data, _freq = self.rebin_freq(data/bdps[None, :, :], n=20)
        #self.plot_wf(_data, self.time[:nblock], freq=_freq/1.e3)
        mask = mask.reshape(self.mask.shape)
        self.mask = mask
        
    
    def plot_spec(self, data=None, plot_indx=False, with_bp=False, suffix=None,
            f_rebin=None):
        
        if data is None:
            data = np.ma.array(self.data)
            data.mask = copy.deepcopy(self.mask)
            #print np.any(data.mask)
        freq = self.freq
        if f_rebin is not None:
            data, freq = self.rebin_freq(data, n=f_rebin)
        if self.cal_on is not None:
            data_cal = data[:, self.cal_on, ...]
            data_cal = data_cal.reshape((-1, ) + data_cal.shape[-2:])
            data = data[:, ~self.cal_on, ...]
        data = data.reshape((-1, ) + data.shape[-2:])
        if with_bp:
            bdps = self.get_bandpass(data, kernel=301)
        #data = data / bdps[None, :, :]
        #poly = self.get_polymean(data, deg=3)
        #data += poly[None, :, :]
        
        spec_mean = np.ma.mean(data, axis=0)
        if self.cal_on is not None:
            spec_cal = np.ma.mean(data_cal, axis=0)
        #spec_mean = np.ma.median(data, axis=0)

        fig  = plt.figure(figsize=(10, 6))
        axhh = fig.add_axes([0.10, 0.52, 0.85, 0.40])
        axvv = fig.add_axes([0.10, 0.10, 0.85, 0.40])
        
        if plot_indx:
            xx = np.arange(freq.shape[0])
        else:
            xx = freq * 1.e-3
        
        axhh.plot(xx, spec_mean[0, ], '-', linewidth=0.8)
        axvv.plot(xx, spec_mean[1, ], '-', linewidth=0.8)
        if self.cal_on is not None:
            axhh.plot(xx, spec_cal[0, ], '-', linewidth=0.8)
            axvv.plot(xx, spec_cal[1, ], '-', linewidth=0.8)
        
        
        #axhh.plot(self.freq * 1e-3, poly[0, ], 'r-', linewidth=1)
        #axvv.plot(self.freq * 1e-3, poly[1, ], 'r-', linewidth=1)
        
        if with_bp:
            axhh.plot(xx * 1e-3, bdps[0, ], 'g-', linewidth=0.8)
            axvv.plot(xx * 1e-3, bdps[1, ], 'g-', linewidth=0.8)
        
        axhh.set_xticklabels([])   
        axhh.set_ylabel('V [XX Polarization]')
        axhh.minorticks_on()
        axhh.tick_params(length=4, width=0.8, direction='in')
        axhh.tick_params(which='minor', length=2, width=0.8, direction='in')
        
        if plot_indx:
            axvv.set_xlabel('Frequency \#')
        else:
            axvv.set_xlabel(r'Frequency [GHz]')
        axvv.set_ylabel('V [YY Polarization]')
        axvv.minorticks_on()
        axvv.tick_params(length=4, width=0.8, direction='in')
        axvv.tick_params(which='minor', length=2, width=0.8, direction='in')

        if suffix is not None:
            fig.savefig('./png/spec_%s.png'%suffix)

        #plt.show()

    def flag_cal(self, p=8, l=1, d=0):
        
        cal_on = np.zeros(self.data.shape[1]).astype('bool')
        cal_on = cal_on.reshape(-1, p)
        cal_on[:, d:l] = True
        cal_on = cal_on.flatten()
        self.cal_on = cal_on
        self.cal_off = np.roll(cal_on, 1)
        #self.mask[:, cal_on,...] = True
        
    def plot_timestream(self, data = None, suffix=None, flag_cal_on=True ):
        
        if data is None:
            data = np.ma.array(self.data)
            data.mask = copy.deepcopy(self.mask)
        if flag_cal_on and self.cal_on is not None:
            data.mask[:, self.cal_on, ...] = True
        #data = self.data
        data = data.reshape((-1, ) + data.shape[-2:])
        time_mean = np.ma.mean(data, axis=-1)
        print time_mean.max()
        print time_mean.min()
        
        
        fig  = plt.figure(figsize=(10, 6))
        axhh = fig.add_axes([0.10, 0.52, 0.85, 0.40])
        axvv = fig.add_axes([0.10, 0.10, 0.85, 0.40])
        
        time = self.time + self.date_obs
        x_axis = [ datetime.fromtimestamp(s.unix) for s in time.flatten()]
        x_label = '%s' % x_axis[0].date()
        x_axis = mdates.date2num(x_axis)
        print x_axis[0], x_axis[-1]
        
        #x_axis = np.arange(time.shape[1])
        
        
        axhh.plot(x_axis, time_mean[:, 0])
        axvv.plot(x_axis, time_mean[:, 1])
        
        date_format = mdates.DateFormatter('%H:%M:%S')
        axhh.xaxis.set_major_formatter(date_format)
        axhh.set_xticklabels([])
        axhh.minorticks_on()
        axhh.tick_params(length=4, width=0.8, direction='in')
        axhh.tick_params(which='minor', length=2, width=0.8, direction='in')
        
        axvv.xaxis.set_major_formatter(date_format)
        axvv.set_xlabel(x_label)
        axvv.minorticks_on()
        axvv.tick_params(length=4, width=0.8, direction='in')
        axvv.tick_params(which='minor', length=2, width=0.8, direction='in')
        
        fig.autofmt_xdate()
        if suffix is not None:
            fig.savefig('./png/timestream_%s.png'%suffix)
        #plt.show()
    

    def cal_to_cal(self, data):

        print data.shape
        cal_on = data[:, self.cal_on, ...]
        cal_off = data[:, self.cal_off, ...]

        cal = np.mean(cal_on - cal_off, axis=1)
        data = data / cal[:, None, ...]

        return data
        
        
    def plot_wf(self, data=None, time=None, freq=None, vmin=None, vmax=None,
                      xmin=None, xmax=None, ymin=None, ymax=None,
                      suffix=None, f_rebin=None, sub_mean=True, 
                      cal_to_cal=True):
        
        fig  = plt.figure(figsize=(10, 6))
        axhh = fig.add_axes([0.10, 0.52, 0.75, 0.40])
        axvv = fig.add_axes([0.10, 0.10, 0.75, 0.40])
        cax  = fig.add_axes([0.86, 0.20, 0.02, 0.60])

        #im = axhh.pcolormesh(x_axis, y_axis, vis1[:,:,0].T, vmax=vmax, vmin=vmin)
        #im = axvv.pcolormesh(x_axis, y_axis, vis1[:,:,1].T, vmax=vmax, vmin=vmin)
        
        
        
        
        if data is None:
            data = np.ma.array(self.data)
            data.mask = copy.deepcopy(self.mask)
        if time is None:
            time = self.time
        if freq is None:
            freq = self.freq / 1.e3
        if f_rebin is not None:
            data, freq = self.rebin_freq(data, n=f_rebin)
            freq = freq / 1.e3
        if cal_to_cal:
            print 'calibrated to the mean noise cal'
            data = self.cal_to_cal(data)
            
        y_axis = freq
        #print y_axis
                       
        
        time = (time + self.date_obs).flatten()
        x_axis = [ datetime.fromtimestamp(s.unix) for s in time.flatten()]
        x_label = '%s' % x_axis[0].date()
        x_axis = mdates.date2num(x_axis)
        #print x_axis
        #print x_label
        
        if self.cal_on is not None:
            data.mask[:, self.cal_on, ...] = True


        data = data.reshape((-1, ) + data.shape[-2:])
        data = data[:, :2, :]

        if sub_mean:
            data -= np.ma.mean(data, axis=(0, -1))[None, :, None]
        
        if vmin is None:
            _m = np.ma.mean(data)
            _v = np.ma.std(data)
            vmin = _m - 3. * _v
            vmax = _m + 3. * _v
            print vmin, vmax
                
        im = axhh.pcolormesh(x_axis, y_axis, data[:, 0, :].T, vmin=vmin, vmax=vmax)
        im = axvv.pcolormesh(x_axis, y_axis, data[:, 1, :].T, vmin=vmin, vmax=vmax)
            
        fig.colorbar(im, cax=cax, ax=axvv)
        
        date_format = mdates.DateFormatter('%H:%M:%S')
        axhh.xaxis.set_major_formatter(date_format)
        axhh.set_xticklabels([])
        axhh.set_ylabel(r'${\rm frequency\, [GHz]\, HH}$')
        if xmin is None: xmin = x_axis[0]
        if xmax is None: xmax = x_axis[-1]
        if ymin is None: ymin = y_axis[0]
        if ymax is None: ymax = y_axis[-1]
        axhh.set_xlim(xmin=xmin, xmax=xmax)
        axhh.set_ylim(ymin=ymin, ymax=ymax)
        axhh.minorticks_on()
        axhh.tick_params(length=4, width=0.8, direction='in')
        axhh.tick_params(which='minor', length=2, width=0.8, direction='in')
        
        axvv.xaxis.set_major_formatter(date_format)
        axvv.set_xlabel(x_label)
        axvv.set_ylabel(r'${\rm frequency\, [GHz]\, VV}$')
        axvv.set_xlim(xmin=xmin, xmax=xmax)
        axvv.set_ylim(ymin=ymin, ymax=ymax)
        axvv.minorticks_on()
        axvv.tick_params(length=4, width=0.8, direction='in')
        axvv.tick_params(which='minor', length=2, width=0.8, direction='in')


        fig.autofmt_xdate()

        #cax.set_ylabel(r'${\rm V}/{\rm V}_{\rm time median}$')
        #cax.set_ylabel(r'${\rm V}/{\rm V}_{\rm noise\, cal}$')
        cax.set_ylabel('V')
        if suffix is not None:
            fig.savefig('./png/wf_%s.png'%suffix)
        #plt.show() 
        
