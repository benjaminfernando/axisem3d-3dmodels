#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from obspy.core import Trace, Stats, UTCDateTime
import matplotlib.pyplot as plt

#you can obviously have as many traces plotted as you want. I've used four, for clarity - obviously you need to alter the paths. 
#Of course, there are much neater ways to write these scripts to loop over the different files, but I think these are the easiest to follow. 

file_blank =  '/Users/benfernando/Desktop/Leeds/data/synthetic_data/D1/copy_off/global_seismic_network_GSN/'
file_topo_only = '/Users/benfernando/Desktop/Leeds/data/synthetic_data/D2/copy_off/global_seismic_network_GSN/'
file_moho_only = '/Users/benfernando/Desktop/Leeds/data/synthetic_data/D5/copy_off/global_seismic_network_GSN/'
file_topo_moho = '/Users/benfernando/Desktop/Leeds/data/synthetic_data/D6/copy_off/global_seismic_network_GSN/'

def load_data(station,component, file):
    
    axisem3d_data = np.loadtxt(file + str(station+'.ascii'))
    times_a = np.loadtxt(file+'data_time.ascii')
    axisem3d_stats = Stats()
    axisem3d_stats.delta = (times_a[1] - times_a[0])
    axisem3d_stats.starttime = UTCDateTime(times_a[0])
    axisem3d_stats.npts = times_a.size
    component_data=np.hsplit(axisem3d_data,3)
    
    if component =='z':
        c_a = np.concatenate(component_data[2])
    elif component =='t':
        c_a = np.concatenate(component_data[1])
    elif component =='r':
        c_a = np.concatenate(component_data[0])
    else:
        raise Exception('No component with this name')
    
    axisem3d_trace = Trace(c_a, header=axisem3d_stats)
    return axisem3d_trace

def comparision_plot(station, component,period):
    
    plt.figure()
    
    trace1 = load_data(station, component, file_blank)
    trace2 = load_data(station, component, file_topo_only)
    trace3 = load_data(station, component, file_moho_only)
    trace4 = load_data(station, component, file_topo_moho)

    
    trace1 = filter_data(trace1,period)
    trace2 = filter_data(trace2,period)
    trace3 = filter_data(trace3,period)
    trace4 = filter_data(trace4,period)

    
    plt.plot(trace1.times(), trace1.data*1e3,label='Without Topography or Moho', lw=.5, color = 'b')
    plt.plot(trace2.times(), trace2.data*1e3,label='With Topography only', lw=.5, color = 'r')
    plt.plot(trace3.times(), trace3.data*1e3,label='With Moho only', lw=.5, color = 'g')
    plt.plot(trace4.times(), trace4.data*1e3,label='With Topography and Moho', lw=.5, color = 'k')

    plt.ylabel('Displacement amplitude (millimetres)')
    plt.xlabel('Time since event (seconds)')
    plt.title(str(station) + " " + str(period) + ' seconds')
    plt.legend()
    plt.show()
    
def filter_data(trace, shortest_period):

    trace = trace.filter('bandpass', freqmax=1.0/shortest_period, freqmin=1.0/100)
    
    return trace