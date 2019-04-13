#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:11:32 2019

@author: jon
"""

from astropy.io import fits
import numpy as np
from msinput import MS_jon

def convert_simulation_point():
    ms = MS_jon()
    ms.read_ms("./ms/simulation_point/simulation.ms", data_column="CORRECTED_DATA")
    
    freq  = ms.freq_array[0]
    baselines = ms.uvw_array.shape[0] // ms.Ntimes
    
    uvw = np.zeros((baselines, ms.Ntimes, 3), dtype=np.double)
    visibilities = np.zeros((baselines, ms.Ntimes, freq.size, 4, 2), dtype=np.double)
    start=0
    end=baselines
    for i in range(0, ms.Ntimes):
        uvw[:, i] = ms.uvw_array[start:end]
        visibilities[:, i, :, :, 0] = np.real(ms.data_array[start:end, 0,])
        visibilities[:, i, :, :, 1] = np.imag(ms.data_array[start:end, 0,])
        start += baselines
        end += baselines
    
    hdulFreq = fits.PrimaryHDU(freq)
    hdulFreq.writeto("fits/simulation_point/freq.fits")
    hdulUVW = fits.PrimaryHDU(uvw)
    hdulUVW.writeto("fits/simulation_point/uvw.fits")
    hdulVis = fits.PrimaryHDU(visibilities)
    hdulVis.writeto("fits/simulation_point/vis.fits")


ms = MS_jon()
ms.read_ms("./ms/meerkat_tiny/splitted.ms", data_column="DATA")
freq  = ms.freq_array[0]
hdulFreq = fits.PrimaryHDU(freq)
hdulFreq.writeto("fits/meerkat_tiny/freq.fits")

vis_block = ms.data_array
vis_block[ms.flag_array] = 0.0

non_autocor = ms.ant_1_array != ms.ant_2_array
vis_block = vis_block[non_autocor]
uvw_block = ms.uvw_array[non_autocor]
flags_block = ms.flag_array[non_autocor]

baselines = uvw_block.shape[0] // ms.Ntimes
splits = 8
split_baselines= baselines // splits
for j in range(0, splits):
    flags = np.zeros((split_baselines, ms.Ntimes, freq.size, 4), dtype=np.double)
    uvw = np.zeros((split_baselines, ms.Ntimes, 3), dtype=np.double)
    visibilities = np.zeros((split_baselines, ms.Ntimes, freq.size, 4, 2), dtype=np.double)
    start=split_baselines*j
    end=split_baselines*(j+1)
    for i in range(0, ms.Ntimes):
        uvw[:, i] = uvw_block[start:end]
        visibilities[:, i, :, :, 0] = np.real(vis_block[start:end, 0,])
        visibilities[:, i, :, :, 1] = np.imag(vis_block[start:end, 0,])
        flags[:, i, :, :] =  flags_block[start:end, 0]
        start += baselines
        end += baselines
    hdulflags = fits.PrimaryHDU(flags)
    hdulflags.writeto("fits/meerkat_tiny/flags"+str(j)+".fits")
    hdulUVW = fits.PrimaryHDU(uvw)
    hdulUVW.writeto("fits/meerkat_tiny/uvw"+str(j)+".fits")
    hdulVis = fits.PrimaryHDU(visibilities)
    hdulVis.writeto("fits/meerkat_tiny/vis"+str(j)+".fits")
    