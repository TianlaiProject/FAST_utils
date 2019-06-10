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

def est_tcorr_psd1d_lombscargle(data, ax, n_bins=None, inttime=None,
        f_min=None, f_max=None):

    fft_len = data.shape[0]

    if inttime is None:
        inttime = ax[1] - ax[0]

    print
    print 'int time', inttime
    print


    #windowf = np.blackman(fft_len)

    if n_bins is None: n_bins = 30
    #if f_min  is None: f_min = 1. / ( ax.max()     - ax.min() )
    #if f_max  is None: f_max = 1. / ( np.abs(ax[1] - ax[0]) ) / 2.
    if f_min  is None: f_min = 1. / ( inttime * fft_len )
    if f_max  is None: f_max = 1. / inttime / 2.

    freq_bins_c = np.logspace(np.log10(f_min), np.log10(f_max), n_bins)
    #freq_bins_c = np.linspace(f_min, f_max, n_bins)
    #freq_bins_d = freq_bins_c[1] / freq_bins_c[0]
    #freq_bins_e = freq_bins_c / (freq_bins_d ** 0.5)
    #freq_bins_e = np.append(freq_bins_e, freq_bins_e[-1] * freq_bins_d)

    #freqs = np.linspace(f_min, f_max, ax.shape[0])
    #freqs = np.linspace(f_min, f_max, 1024)

    n_freq, n_pol = data.shape[1:]

    power = np.zeros([n_bins, n_freq, n_pol])

    #norm = np.histogram(freqs, bins=freq_bins_e)[0] * 1.
    #norm[norm==0] = np.inf
    for i in range(n_freq):
        for j in range(n_pol):
            y = data[:, i, j]
            y = y - np.ma.mean(y)
            #y = y * windowf
            power[:, i, j] = lombscargle(ax, y, 2. * np.pi * freq_bins_c,
                                         normalize=False)
            #hist   = np.histogram(freqs, bins=freq_bins_e, weights=_p)[0]
            #power[:, i, j] = hist / norm
    #power = np.sqrt(4. * power / float(ax.shape[0]) / np.std(y) ** 2.)
    power *= inttime

    power[freq_bins_c <= 1. / ( inttime * fft_len ) ] = 0.
    power[freq_bins_c >= 1. / inttime / 2.] = 0.

    return power, freq_bins_c

def est_fnps(data, time, df=100, output_path = None):
    
    time = time.to(u.s).value
    
    mask = data.mask
    good_freq = ~np.all(mask, axis=(0, 1))
    
    data = data[:, :, good_freq]
    data_shp = data.shape
    print data_shp
    nfreq = data.shape[-1]
    nfreq = int(nfreq / df) * df
    data = data[:, :, :nfreq]
    data = data.reshape(data_shp[:-1] + (-1, df))
    data = np.ma.mean(data, axis=-1)
    print data.shape
    
    power, bins_cent = est_tcorr_psd1d_lombscargle(data, time)
    
    if output_path is not None:
        with h5py.File(output_path, 'w') as fh:
            fh['tcorr_ps'] = power
            fh['tcorr_bc'] = bins_cent
    else:
        return power, bins_cent

def plot_result_one(data_base, data_file, label=None,
                c='r', fmt='o-', mfc='none', axes=None, shift=1.):
    if axes is None:
        fig = plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.05)
        axhh = fig.add_subplot(gs[0, 0])
        axvv = fig.add_subplot(gs[0, 1])
        axes = [axhh, axvv]
    else:
        axhh, axvv = axes
    pols = ['HH', 'VV']
    
    legend_list = []
    
    with h5py.File(data_base + data_file, 'r') as fh:

        tcorr_ps = fh['tcorr_ps'][:]
        tcorr_bc = fh['tcorr_bc'][:]
        
        shift = (((tcorr_bc[1] / tcorr_bc[0]) ** 0.5) ** 0.3)**shift
        
        n_ps = float(tcorr_ps.shape[-1])
        mean = np.mean(tcorr_ps, axis=-1)
        #mean /= mean.max()
        erro = np.std(tcorr_ps, axis=-1)
        erro /= np.sqrt(n_ps)
        
        upper = mean + erro
        lower = mean - erro

        errors = erro
        
        pp = mean[:,0] > 0
        #axhh.errorbar(tcorr_bc[pp] * shift, mean[pp, 0], errors[pp, 0], fmt=fmt[0],
        #        c=c, marker=fmt[0], mfc=mfc, mec=c, ms=5, mew=1, lw=1,
        #        ecolor=c, elinewidth=1.5, capsize=3, capthick=1.5)
        axhh.plot(tcorr_bc[pp] * shift, erro[pp, 0], color=c, linewidth=1.5,
                  drawstyle='steps-mid')
        
        pp = mean[:,1] > 0
        #axvv.errorbar(tcorr_bc[pp] * shift, mean[pp, 1], errors[pp, 1], fmt=fmt[0],
        #        c=c, marker=fmt[0], mfc=mfc, mec=c, ms=5, mew=1, lw=1,
        #        ecolor=c, elinewidth=1.5, capsize=3, capthick=1.5)
        axvv.plot(tcorr_bc[pp] * shift, erro[pp, 1], color=c, linewidth=1.5,
                  drawstyle='steps-mid')
        
    if (label is not None) and (label is not ''):
        legend_list.append(mpatches.Patch(color=c, label=label))
    xmin = tcorr_bc.min() / 1.5
    xmax = tcorr_bc.max() * 1.8
    
    return xmin, xmax, legend_list

def plot_results_group(data_base, data_group_list, label_list=None, vmin=None, vmax=None,
                       tcorr=True, output_path=None, title=''):
    
    fig = plt.figure(figsize=(8, 5))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.02, top=0.93)
    axhh = fig.add_subplot(gs[0, 0])
    axvv = fig.add_subplot(gs[0, 1])
    axes = [axhh, axvv]
    
    legend_list = []
    xmin = 1.e20
    xmax = -1.e20
    
    shift = np.arange(len(data_group_list)) - (len(data_group_list) - 1) * 0.5
    
    for ii, data_group in enumerate(data_group_list):
        
        c = _c_list[ii]
        mfc = c

        for jj, data_file in enumerate(data_group):
            if mfc is not 'w': mfc = 'w'
            else: mfc = c
                
            if label_list is not None: label = label_list[ii][jj]
                
            fmt = _l_list[jj]
            o = plot_result_one(data_base, data_file, label=label,
                c=c, fmt=fmt, mfc=mfc, axes=axes, shift=shift[ii])
            if xmin > o[0]: xmin = o[0]
            if xmax < o[1]: xmax = o[1]
            legend_list += o[2]
    
    locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                                      numticks=100)
    if tcorr:
        axhh.set_xlabel(r'$f\,[{\rm Hz}]$')
        #axhh.set_ylabel(r'${\rm PSD}(f)$')
        axhh.set_ylabel(r'$S^{t}(f)\,[{\rm s}^{-1}]$')
    else:
        axhh.set_xlabel(r'$\omega\,[1/{\rm MHz}]$')
        #axhh.set_ylabel(r'${\rm PSD}(\omega)$')
        axhh.set_ylabel(r'$S^{\nu}(\omega)\,[{\rm MHz}^{-1}]$')
        
    #axhh.set_title(r'${\rm HH\,\, Pol.}$')
    if label is not None:
        _l = axhh.legend(handles=legend_list, frameon=False,
                    #title=r'${\rm HH\,\, Pol.}$', 
                    title = title + 'HH Polarization',
                    markerfirst=False)
        _l._legend_box.align = "right"
        _l = axvv.legend(frameon=False,
                    #title=r'${\rm VV\,\, Pol.}$', 
                    title= title + 'VV Polarization',
                    markerfirst=False)
        _l._legend_box.align = "right"
    axhh.set_xlim(xmin=xmin, xmax=xmax)
    axhh.set_ylim(ymin=vmin, ymax=vmax)
    axhh.loglog()

    axhh.tick_params(length=4, width=1, direction='in')
    axhh.tick_params(which='minor', length=2, width=1, direction='in')
    axhh.xaxis.set_major_locator(locmaj)
    axhh.xaxis.set_minor_locator(locmin)
    axhh.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    if tcorr:
        axvv.set_xlabel(r'$f\,[{\rm Hz}]$')
        fig.suptitle('Temporal Power Spectral Density')
    else:
        axvv.set_xlabel(r'$\omega\,[1/{\rm MHz}]$')
        fig.suptitle('Spectroscopic Power Spectral Density')
    #axvv.set_title(r'${\rm VV\,\, Pol.}$')
    axvv.set_xlim(xmin=xmin, xmax=xmax)
    axvv.set_ylim(ymin=vmin, ymax=vmax)
    axvv.loglog()
    #axvv.semilogx()
    axvv.set_yticklabels([])
    #axvv.minorticks_on()
    axvv.tick_params(length=4, width=1, direction='in')
    axvv.tick_params(which='minor', length=2, width=1, direction='in')
    axvv.xaxis.set_major_locator(locmaj)
    axvv.xaxis.set_minor_locator(locmin)
    axvv.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    if output_path is not None:
        fig.savefig(output_path + '.eps', formate='eps')
        
    plt.show()
    fig.clf()
    plt.close()
    
