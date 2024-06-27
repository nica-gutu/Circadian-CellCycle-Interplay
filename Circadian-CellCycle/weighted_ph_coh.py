#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:59:24 2023

@author: nicagutu
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em
from pyboat import plotting as pl 
import math

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def func_weight(k, dist):
    return np.exp(-k * dist)

def normalize(signal):
    return (signal-min(signal))/(max(signal)-min(signal))

dt = 0.5 # the sampling interval, 0.5hours
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

path = 'Data/'
output = 'Output/'

cond = '100'
cell_annot = pd.read_excel(path+r'cell_annotations_'+str(cond)+'.xlsx', header=None)
list_pos = cell_annot.index[cell_annot[1].diff() != 0].tolist()
list_pos[0] = 0
list_pos.append(cell_annot.index.values[-1]+1)
print(list_pos)

factor = 6.5*2/20 #from pixels to uM

k_vect = [0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 1] ###exponential decay of the weight function

ph_coh_position = pd.DataFrame()
for j in range(len(list_pos)-1):
    start = list_pos[j]
    end = list_pos[j+1]
    print(start, end)
        
    file = 'Circadian_traces_'+str(cond)
    data = pd.read_excel(path+file+r'.xlsx',header=None)
    data = data.T
        
    file1 = 'Centroid_row_'+str(cond)
    file2 = 'Centroid_col_'+str(cond)
    data_row = pd.read_excel(path+file1+r'.xlsx', header=None)
    data_col = pd.read_excel(path+file2+r'.xlsx', header=None)
    data_row = data_row.T
    data_col = data_col.T
    
    data = data.iloc[:, start:end]

    time = data.index.values*0.5
    
    data_row = data_row.iloc[:, start:end]*factor
    data_col = data_col.iloc[:, start:end]*factor
    
    phases = pd.DataFrame()
    ridge_results = {}
    filtered = []
    for col in data:
        signal = data[col]
        # if len((signal.replace(-1,np.nan)).dropna())>150:
        detrended = wAn.sinc_detrend(signal, T_c=50)
        norm = wAn.normalize_amplitude(detrended, window_size=50) 
        wAn.compute_spectrum(detrended, do_plot=False)
        rd = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)    
        ridge_results[col] = rd
        phases[col] = rd['phase']
        
    # powers_series = em.average_power_distribution(ridge_results.values(), signal_ids = ridge_results.keys())
    # high_power_ids = powers_series[powers_series > 0].index
    # high_power_ridge_results = [ridge_results[i] for i in high_power_ids]
    # res = em.get_ensemble_dynamics(high_power_ridge_results)
        
    # for col in data:
    #     signal = normalize(data[col])
    #     detrended = wAn.sinc_detrend(signal, T_c=50)
    #     norm = wAn.normalize_amplitude(detrended, window_size=50) 
    #     smoothed = wAn.sinc_smooth(norm, T_c=20)
    #     plt.plot(time, normalize(smoothed), alpha=0.5)   
    # plt.xlabel('Time [hours]')
    # plt.xticks([0,24,48,72,96,120])
    # plt.ylabel('Signal [a.u.]')
    # plt.savefig(output+'Signals_position'+str(j)+'.svg')
    # plt.show()
        
    plt.figure(figsize=(12,10))
    # plt.plot(time, res[3]['R'], linewidth=5, label='Standard')    
    for k in k_vect:    
        ph_coh_df_loc = pd.DataFrame()
        for col in data:
            cell = data[col]
            col_ref, row_ref = data_col[col], data_row[col]

            # data_new = data.drop(col, axis=1)               
            R_all = np.zeros(len(time))
            for ii in range(len(time)):
                sum_phases_inst = 0
                count = 0
                for column in data:
                    col_a, row_a = data_col[column][ii], data_row[column][ii]
                    dist = np.sqrt((col_a-col_ref[ii])**2+(row_a-row_ref[ii])**2)
                    phase_eff = (phases[column][ii])
                    sum_phases_inst = sum_phases_inst+np.exp(phase_eff*1j)*func_weight(k, dist) #-phases[cell][ii]????
                    count += np.median(func_weight(k, dist))
                R_all[ii] = np.abs(sum_phases_inst/count)
            R_all[len(R_all)-1] = R_all[len(R_all)-2]
            ph_coh_df_loc[col] = R_all
        
        if k==0:
            ph_coh_position[j] = ph_coh_df_loc.mean(axis=1)

        plt.plot(time, ph_coh_df_loc.mean(axis=1), linewidth=5, label=str(k))
    plt.xlabel('Time [hours]')
    plt.xticks([0,24,48,72,96,120])
    plt.ylabel('Phase coherence')
    plt.legend(loc='best')
    plt.savefig(output+'Phase_coh_decay_position'+str(j)+'.svg')
    plt.show()
            

file = 'Circadian_traces_'+str(cond)
data = pd.read_excel(path+file+r'.xlsx',header=None)
data = data.T
for col in data:
    signal = data[col]
    # if len((signal.replace(-1,np.nan)).dropna())>150:
    detrended = wAn.sinc_detrend(signal, T_c=50)
    norm = wAn.normalize_amplitude(detrended, window_size=50) 
    wAn.compute_spectrum(detrended, do_plot=False)
    rd = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)    
    ridge_results[col] = rd
    phases[col] = rd['phase']

powers_series = em.average_power_distribution(ridge_results.values(), signal_ids = ridge_results.keys())
high_power_ids = powers_series[powers_series > 0].index
high_power_ridge_results = [ridge_results[i] for i in high_power_ids]
res = em.get_ensemble_dynamics(high_power_ridge_results)

plt.figure(figsize=(12,10))
plt.plot(time, ph_coh_position.mean(axis=1), linewidth=5, label='Mean positions')
plt.plot(time, res[3]['R'], linewidth=5, label='Standard')
plt.xlabel('Time [hours]')
plt.xticks([0,24,48,72,96,120])
plt.ylabel('Phase coherence')
plt.legend(loc='best')
plt.savefig(output+'Standard_vs_mean_positions.svg')
plt.show()

dimx = 1022*factor
dimy = 1024*factor

total_dist = np.arange(0, np.sqrt(dimx**2+dimy**2), step=1)
plt.figure(figsize=(12,10))
for k in k_vect:
    plt.plot(total_dist, func_weight(k, total_dist), linewidth=5, label=str(k))
plt.xlabel('Distance [uM]')
plt.ylabel('Weight function')
legend = plt.legend(loc='best')
legend.set_title('Decay')
plt.savefig(output+'Weight_function_decay_one_position.svg')
plt.show()    
    




