#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 11:36:02 2019

Note:

    
_______________________________________________________________________________

author: 
    @71unxv
    hibatullah.geophysics@gmail.com
    github.com/71unxv
    linkedin.com/in/irsyadhbtlh/
    
Cheers!
"""
###############################################################################
# Import Library
import numpy as np
import bruges.filters.wavelets as wv
from psopy import _minimize_pso

###############################################################################
# Function
def objWavelet(phase):
    global freq,RC,seis,w
    
    wavelet = change_phase(w,phase,dt, degrees=False)
    syn_seis = np.convolve(RC,wavelet, mode='same')
    koef_corr = np.corrcoef(seis,syn_seis)
    
    
    return koef_corr

def change_phase(w,phase, dt, degrees=False):
#    seismic = np.load('wavelet-est-SEG/seismic.npy')

    spec = np.fft.rfft(w, len(w))

    comp = spec * np.exp(1j*phase)
    w_new= np.fft.ifft(comp, n=len(w))
    return w_new
###############################################################################
# Run code! Run!!
global freq,RC,seis,w

RC = np.load("wavelet-est-SEG/rpp.npy")
seis = np.load("wavelet-est-SEG/seismic.npy")
dt = .002  # sample rate in seconds
t = np.arange(seis.size) * dt

###############################################################################
# create wavelet

freqs = np.array([5, 8, 130, 160])
w = wv.ormsby(0.5, dt, freqs)

x0 = np.deg2rad(np.arange(0,90,len(w)))

res = _minimize_pso(objWavelet, x0)


