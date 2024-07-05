#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:07:06 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def rescale(signal):
    minimum = np.min(signal)
    maximum = np.max(signal)
    return (signal-minimum)/(maximum-minimum)

def linear_slope(vector):
    x_values = np.arange(len(vector))
    y_values = np.array(vector)
    slope, intercept = np.polyfit(x_values, y_values, 1)
    return slope

path = 'Raw_data/'
output = 'Circ_CellCycle/'
dose = ['untreated', '5uM', '10uM'] #'untreated',, '1.25uM','2.5uM',
density = ['high'] #,'medium','low'
channel1 = 'circadian'
channel2 = 'cell_cycle'

dt = 0.5 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

for i in density:
    for j in dose:

        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Circadian data
        data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
        #Cell-cycle data
        data2 = pd.read_csv(path+condition+channel2+r'_filtered.csv', index_col=0)        
        
        count = 0
        count2 = 0
        max_count = 12
        
        pd_circ_signals = {}
        pd_circ_time = {}
        
        pd_cc_signals = {}
        pd_cc_time = {}
        
        for column in data1:
            norm_signal = data1[column].dropna()
            
            signal = data2[column].dropna()
            norm_signal2 = wAn.normalize_amplitude(signal, window_size=50)
            
            peaks1,_ = find_peaks(norm_signal, distance=40)
            peaks2,_ = find_peaks(norm_signal2, distance=40)

            wAn.compute_spectrum(norm_signal, do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            power = np.mean(rd['power'])

            wAn.compute_spectrum(norm_signal2, do_plot=False)
            rd2 = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            power2 = np.mean(rd['power'])
            
            if 0<=count<max_count and len(peaks1)>=4 and power>20:
                pd_circ_time[column] = norm_signal.index.values*0.5
                pd_circ_signals[column] = norm_signal
                count += 1
                
            if 0<=count2<max_count and len(peaks2)>=4 and power>10:
                pd_cc_time[column] = norm_signal2.index.values*0.5
                pd_cc_signals[column] = rescale(norm_signal2)
                count2 += 1                 
        
        fig, axes = plt.subplots(1, max_count, figsize=(5 * count, 6), sharey=True)
        for ax, column in zip(axes, pd_circ_signals.keys()):
            ax.plot(pd_circ_time[column], pd_circ_signals[column], linewidth=5)
            ax.set_xlabel('Time [hours]')
        
        axes[0].set_ylabel('Detrended Circadian signal')          
        plt.tight_layout()
        plt.savefig(output+'Circadian_signals_'+str(j)+'.svg')
        plt.show()
        
        fig, axes = plt.subplots(1, max_count, figsize=(5 * count2, 6), sharey=True)
        for ax, column in zip(axes, pd_cc_signals.keys()):
            ax.plot(pd_cc_time[column], pd_cc_signals[column], linewidth=5)
            ax.set_xlabel('Time [hours]')
        
        axes[0].set_ylabel('Detrended CC signal')          
        plt.tight_layout()
        plt.savefig(output+'CC_signals_'+str(j)+'.svg')
        plt.show()
                
            
            
            # wAn.compute_spectrum(norm_signal, do_plot=False)
            # rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            # diff_rel = (rd['phase']-rd2['phase'])
            # diff = np.arctan2(np.sin(diff_rel), np.cos(diff_rel))

            # rand = np.random.uniform(0,1)

            # if rand>0.5 and count<20 and (np.std(diff)<=0.5) and (len(peaks1) and len(peaks2))>=4:
            #     time = norm_signal.index.values
            #     plt.figure(figsize=(10,8))
            #     plt.plot(time, rescale(norm_signal), linewidth=5, color='tab:brown', label='circadian clock')
            #     plt.plot(time, rescale(norm_signal2), linewidth=5, color='tab:green', label='cell cycle')
            #     plt.xlabel('Time [hours]')
            #     plt.ylabel('Circadian clock signal [a.u.]')
            #     # plt.savefig(output+'Signal_'+str(count)+'.svg')
            #     plt.show()

            #     # plt.figure(figsize=(10,8))
            #     # plt.plot(time, norm_signal2, linewidth=5, color='tab:green', label='cell cycle')
            #     # plt.xlabel('Time [hours]')
            #     # plt.ylabel('Cell-cycle signal [a.u.]')
            #     # plt.savefig(output+'Cellcycle_signal_'+str(count)+'.svg')
            #     # plt.show()
                
            #     count+=1

            
            