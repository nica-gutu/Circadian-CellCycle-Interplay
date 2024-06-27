#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:00:17 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from scipy.signal import find_peaks

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path = '.../Raw_data/'
output = '...'
dose = ['untreated', '5uM', '10uM']
density = ['high']
channel1 = 'circadian'
channel2 = 'cell_cycle'
colors = ['gold', 'tab:green', 'tab:blue']

dt = 0.5 
lowT = 16
highT = 40
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#%%Changing dose within same density

num_peaks_cond = []
for i in density:
      
    for j in dose:
        count = 0

        condition = str(i)+'_density_'+str(j)+'_'
        
        #Circadian data
        data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
            
        time = data1.index.values*0.5
        
        for column in data1:
            norm_signal = data1[column].dropna()
            
            if len(norm_signal) > 0:
                wAn.compute_spectrum(norm_signal,do_plot=False)
                rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
                rd.set_index(norm_signal.index,inplace=True)
                mean_power = (np.mean(rd['power']))
                
                peaks = find_peaks(norm_signal, height=0, distance=32, prominence=0.1)
                
                if len(peaks[0]) > 1 and mean_power > 10:
                    count += 1
                
        
        num_peaks_cond.append((count/len(data1.columns))*100)
           
                           
    print(num_peaks_cond)
    
    fig = plt.figure(figsize=(10,10))
    barlist = plt.bar(dose, num_peaks_cond)
    for ii in range(len(dose)):
        barlist[ii].set_color(colors[ii])
    plt.ylabel('% oscillating cells')
    plt.ylim([0,100])
    plt.savefig(output+'Number_osc_cells_high_density.svg')
    plt.show()



#%% Changing density        
density = ['high', 'medium', 'low']

num_peaks_cond = []
for jj in density:
    count = 0 

    condition = str(jj)+'_density_'+str(dose[0])+'_'
    print(condition)
    
    #Circadian data
    data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
        

    for column in data1:
        norm_signal = data1[column].dropna()

        if len(norm_signal) > 0:
            wAn.compute_spectrum(norm_signal,do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            rd.set_index(norm_signal.index,inplace=True)
            mean_power = (np.mean(rd['power']))
            
            peaks = find_peaks(norm_signal, height=0, distance=32, prominence=0.1)
            
            if len(peaks[0]) > 1 and mean_power > 10:
                count += 1
            
    
    num_peaks_cond.append((count/len(data1.columns))*100)
       
                       
print(num_peaks_cond)

fig = plt.figure(figsize=(10,10))
barlist = plt.bar(density, num_peaks_cond)
for ii in range(len(density)):
    barlist[ii].set_color(colors[ii])
plt.ylabel('% oscillating cells')
plt.ylim([0,100])
plt.savefig(output+'Number_osc_cells_condition_changing_density.svg')
plt.show()

