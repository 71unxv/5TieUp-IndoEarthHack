#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 13:59:13 2019

Note:
coba plot antara syntethic vs true
    
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
import matplotlib.pyplot as plt
import lasio as ls


###############################################################################
# Function
def createRC(Z):
    RC = (Z[1:] - Z[:-1]) / (Z[1:] + Z[:-1])
    RC = np.append([0],RC)
    return RC

def ricker(f, length, dt):
    t = np.linspace(-length / 2, (length-dt) / 2, length / dt)
    y = (1. - 2.*(np.pi**2)*(f**2)*(t**2))*np.exp(-(np.pi**2)*(f**2)*(t**2))
    return y
###############################################################################
# Run code! Run!!


well_raw = ls.read('5tieup/Kapuni-8/Wireline/kapuni-8_final.las')
well_data=well_raw.df()
DEPTH_ALL=np.array(well_data.index)
DENS = np.array(well_data[well_data["DENS"].notnull() & well_data["DTC"].notnull()]["DENS"])
SLOW  = np.array(well_data[well_data["DENS"].notnull() & well_data["DTC"].notnull()]["DTC"])
DEPTH = np.array(well_data[well_data["DENS"].notnull() & well_data["DTC"].notnull()].index)
###########################################################################
# load check shot

chk = np.loadtxt('chk-kapuni-8.txt')
chk_depth =chk[:,0]
chk_time  =chk[:,1] 

#############################################################################
# load seismic traces

seis = np.loadtxt('trace_seismik_kapuni8-1.dat')
seis_time =seis[:,0]/1000
seis_amplitude  =seis[:,1] 

depth_sample = 0.1524
time_sample = 0.002
time_length = 3
ft2m = 0.3048

##############################################################################
# get the T0 value
chk_func = np.poly1d(np.polyfit(chk_depth, chk_time, 2))

T0=chk_func(DEPTH[0])

############################################################################
# calculate Z dan RC

Z = SLOW * DENS
RC = createRC(Z)

######################################################################33
# create synthetic seismic
DT = SLOW * depth_sample *(1/304800)

TIME = 2* np.cumsum(DT)
TIME = TIME + (2*T0)

TIME_resample = np.arange(0, time_length, time_sample) 

Z_resample = np.interp(TIME_resample, TIME, Z)
RC_resample = createRC(Z_resample)


wavelet = ricker(20, 0.512, 0.001)
syn_seismic = np.convolve(RC_resample, wavelet, mode='same')

##############################################################################
# create time by sonic log
well_tie_func = np.poly1d(np.polyfit(DEPTH, TIME, 5))

TIME_OUT=well_tie_func(DEPTH_ALL)

############################################################################
# plot the well
fig01=plt.figure()

fig01.add_subplot(171)

fig01.add_subplot(171).plot(DENS,-TIME)
fig01.add_subplot(171).set_ylim(-2.1,-1.4)


fig01.add_subplot(172)
fig01.add_subplot(172).plot(SLOW,-TIME)
fig01.add_subplot(172).set_ylim(-2.1,-1.4)

fig01.add_subplot(173)
fig01.add_subplot(173).plot(Z,-TIME)
fig01.add_subplot(173).set_ylim(-2.1,-1.4)

fig01.add_subplot(174)
fig01.add_subplot(174).plot(RC,-TIME)
fig01.add_subplot(174).set_ylim(-2.1,-1.4)

fig01.add_subplot(175)
fig01.add_subplot(175).plot(RC_resample,-TIME_resample)
fig01.add_subplot(175).set_ylim(-2.1,-1.4)

fig01.add_subplot(176)
fig01.add_subplot(176).plot(syn_seismic,-TIME_resample)
fig01.add_subplot(176).set_ylim(-2.1,-1.4)

fig01.add_subplot(177)
fig01.add_subplot(177).plot(seis_amplitude,-seis_time)
fig01.add_subplot(177).set_ylim(-2.1,-1.4)

#correlation_seis=np.corrcoef(syn_seismic,seis_amplitude[:1500])
#fig01.add_subplot(177).set_xlim(-4000,4000)
fig01.show()

#fig02=plt.figure()
#
#fig02.add_subplot(111).plot(DEPTH,TIME,'.')
#fig02.show()

