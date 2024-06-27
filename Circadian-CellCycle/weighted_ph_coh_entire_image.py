#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:10:31 2023

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
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def func_weight(k, dist):
    return np.exp(-k * dist)

def func_fit(x, c, r, x0, a):
    return c/(1+np.exp(r*(x-x0)))+a

def func_fit2(x, a, b):
    return a*x+b

def exp_func(k, dist):
    # k_array = np.array(k, dtype=float)
    # dist_array = np.array(dist, dtype=float)

    return 0.08*np.exp(k*dist)

dt = 0.5 # the sampling interval, 0.5hours
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

path = 'Data1/'
path2 = 'Data2/'
output = 'Output/'

# cond = '100'
# cell_annot = pd.read_excel(path+r'cell_annotations_'+str(cond)+'.xlsx', header=None)
# list_pos = cell_annot.index[cell_annot[1].diff() != 0].tolist()
# list_pos[0] = 0
# list_pos.append(cell_annot.index.values[-1]+1)
# print(list_pos)

k_vect = [0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 1] ###exponential decay of the weight function

cond = '100'

#%%Weight function vs decay

total_dist = np.arange(0, np.sqrt(5034**2+5985**2), step=1)
plt.figure(figsize=(12,10))
for k in k_vect:
    plt.plot(total_dist, func_weight(k, total_dist), linewidth=5, label=str(k))
plt.xlabel('Distance [uM]')
plt.xscale('log')
plt.ylabel('Weight function')
legend = plt.legend(loc='best')
legend.set_title('Decay')
plt.savefig(output+'Weight_function_decay_whole_image.svg')
plt.show()   
    
factor = 6.5*2/20
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

#%%Circadian

file = 'Circadian_traces_'+str(cond)
data = pd.read_excel(path+file+r'.xlsx',header=None)
data = data.T
time = data.index.values*0.5

file1 = 'Absolute_centroid_row_'+str(cond)
file2 = 'Absolute_centroid_col_'+str(cond)
data_row = pd.read_excel(path+file1+r'.xlsx', index_col=0, header=0)
data_col = pd.read_excel(path+file2+r'.xlsx', index_col=0, header=0)

phases = pd.DataFrame()
ridge_results = {}

for col in data:
    signal = data[col]
    detrended = wAn.sinc_detrend(signal, T_c=50)
    norm = wAn.normalize_amplitude(detrended, window_size=50) 
    wAn.compute_spectrum(norm, do_plot=False)
    rd = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)    

    phases[col] = rd['phase'] 
    
powers_series = em.average_power_distribution(ridge_results.values(), signal_ids = ridge_results.keys())
high_power_ids = powers_series[powers_series > 0].index
high_power_ridge_results = [ridge_results[i] for i in high_power_ids]
res = em.get_ensemble_dynamics(high_power_ridge_results)

df_ph_coh_decays = pd.DataFrame()

plt.figure(figsize=(12,10))
plt.plot(time, res[3]['R'], linewidth=5, label='Standard')    
for k in k_vect:    
    ph_coh_df_loc = pd.DataFrame()
    print(k)
    for col in phases:
        cell = data[col]
        col_ref, row_ref = data_col[col], data_row[col]

        R_all = np.zeros(len(time))
        for ii in range(len(time)):
            sum_phases_inst = 0
            count = 0
            for column in phases:
                col_a, row_a = data_col[column][ii], data_row[column][ii]
                dist = np.sqrt((col_a-col_ref[ii])**2+(row_a-row_ref[ii])**2)
                phase_eff = (phases[column][ii])
                sum_phases_inst = sum_phases_inst+np.exp(phase_eff*1j)*func_weight(k, dist) #-phases[cell][ii]????
                count += np.median(func_weight(k, dist))
            R_all[ii] = np.abs(sum_phases_inst/count)
        R_all[len(R_all)-1] = R_all[len(R_all)-2]
        R_all[-1] = R_all[-2]
        ph_coh_df_loc[col] = R_all
    
    df_ph_coh_decays[k] = ph_coh_df_loc.mean(axis=1)
    plt.plot(time, ph_coh_df_loc.mean(axis=1), linewidth=5, label=str(k))
plt.xlabel('Time [hours]')
plt.ylabel('Phase coherence of phase differences') #or only phase coherence for circadian
plt.legend(loc='best')
plt.savefig(output+'Phase_coh_decay_entire_image.svg')
plt.show()
        
df_ph_coh_decays.to_excel(path2+'Weighted_ph_coh_entire_image.xlsx')

#%%Circadian
df_ph_coh_decays = pd.read_excel(path2+'Weighted_ph_coh_entire_image.xlsx', index_col=0, header=0)

slope = []
slope_error = []
rel_final_ph_coh = []
decays = [0, 0.0001, 0.001, 0.005, 0.01, 0.02]

# plt.figure(figsize=(12,10))

for col in decays:
    ph_coh = df_ph_coh_decays[col] 
    time = df_ph_coh_decays.index.values*0.5
    params, covariance = curve_fit(func_fit, time, ph_coh, method='lm',  maxfev=10000)   
    # print(params)
    y_fit = func_fit(time, *params)
    # plt.plot(time, df_ph_coh_decays[col], linewidth=5, label=str(col))
    # plt.plot(time, y_fit, '--', linewidth=5)
    rel_final_ph_coh.append(ph_coh.iloc[-1]/ph_coh[0])
    slope.append(params[1])
    slope_error.append(np.sqrt(np.diag(covariance))[1])

              
# plt.xlabel('Time [hours]')    
# plt.ylabel('Phase coherence') 
# plt.xticks([0,24,48,72,96,120])   
# plt.legend(loc='best')
# plt.savefig(output+'Phase_coh_decay_entire_image_fitting.svg')
# plt.show()

plt.figure(figsize=(12,10))
plt.errorbar(decays, slope, yerr=slope_error, fmt='o', markersize=10)
params, covariance = curve_fit(exp_func, decays, slope, method='lm',  maxfev=10000) 
error = np.sqrt(np.diag(covariance))[0]
# print(params)
y_fit2 = exp_func(params[0], np.linspace(0.0001, 0.02, 100))
plt.plot(np.linspace(0.0001, 0.02, 100), y_fit2, '--', linewidth=5, label=r'$e^{{({:.2f} \pm {:.2f}) \cdot dist}}$'.format(params[0], error))
plt.xscale('log', base=10)
# plt.yscale('log', base=10)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Decay rate')
plt.ylabel('Rate phase decoherence')
plt.legend(loc='best')
plt.savefig(output+'Rate_decoherence_circ_fitting.svg')
plt.show()

plt.figure(figsize=(12,10))
plt.plot(decays, rel_final_ph_coh, 'o', markersize=10)
plt.xscale('log', base=10)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Decay rate')
plt.ylabel('Relative phase coherence (final/initial)')
plt.savefig(output+'Relative_decoherence_circ_fitting.svg')
plt.show()

#%%Phase differences
file = 'Circadian_traces_'+str(cond)
data = pd.read_excel(path+file+r'.xlsx',header=None)
data = data.T
time = data.index.values*0.5

file2 = 'Cell_cycle_traces_'+str(cond)
data2 = pd.read_excel(path+file2+r'.xlsx',header=None)
data2 = data2.T

data = data.iloc[:, 0:len(data2.columns)] #Full without cell cycle

file1 = 'Absolute_centroid_row_'+str(cond)
file2 = 'Absolute_centroid_col_'+str(cond)
data_row = pd.read_excel(path+file1+r'.xlsx', index_col=0, header=0)
data_col = pd.read_excel(path+file2+r'.xlsx', index_col=0, header=0)

data_row = data_row.iloc[:, 0:len(data2.columns)] #Full without cell cycle
data_col = data_col.iloc[:, 0:len(data2.columns)] #Full without cell cycle

phases = pd.DataFrame()
ridge_results = {}

power1 = []
power2 = []
for col in data:
    signal = data[col]
    signal2 = data2[col]
    # if len((signal.replace(-1,np.nan)).dropna())>150:
    detrended = wAn.sinc_detrend(signal, T_c=50)
    norm = wAn.normalize_amplitude(detrended, window_size=50) 
    wAn.compute_spectrum(norm, do_plot=False)
    rd = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)    
    power1.append(np.mean(rd['power']))
    
    wAn.compute_spectrum(signal2, do_plot=False)
    rd2 = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)    
    power2.append(np.mean(rd2['power']))
    
    if np.mean(rd2['power'])>12 and np.mean(rd['power'])>15: #for phase differences
        diff = rd['phase']-rd2['phase'] #which oscillator
        phases[col] = np.arctan2(np.sin(diff), np.cos(diff))

        ridge_results[col] = rd
        ridge_results[col]['phase'] = phases[col]
    
    # phases[col] = rd['phase'] #rd only for circadian
    
powers_series = em.average_power_distribution(ridge_results.values(), signal_ids = ridge_results.keys())
high_power_ids = powers_series[powers_series > 0].index
high_power_ridge_results = [ridge_results[i] for i in high_power_ids]
res = em.get_ensemble_dynamics(high_power_ridge_results)

# plt.hist(power1)
# plt.show()
# plt.hist(power2)
# plt.show()    

df_ph_coh_decays = pd.DataFrame()

plt.figure(figsize=(12,10))
plt.plot(time, res[3]['R'], linewidth=5, label='Standard')    
for k in k_vect:    
    ph_coh_df_loc = pd.DataFrame()
    print(k)
    for col in phases:
        cell = data[col]
        col_ref, row_ref = data_col[col], data_row[col]

        # data_new = data.drop(col, axis=1)               
        R_all = np.zeros(len(time))
        for ii in range(len(time)):
            sum_phases_inst = 0
            count = 0
            for column in phases:
                col_a, row_a = data_col[column][ii], data_row[column][ii]
                dist = np.sqrt((col_a-col_ref[ii])**2+(row_a-row_ref[ii])**2)
                phase_eff = (phases[column][ii])
                sum_phases_inst = sum_phases_inst+np.exp(phase_eff*1j)*func_weight(k, dist) #-phases[cell][ii]????
                count += np.median(func_weight(k, dist))
            R_all[ii] = np.abs(sum_phases_inst/count)
        R_all[len(R_all)-1] = R_all[len(R_all)-2]
        R_all[-1] = R_all[-2]
        ph_coh_df_loc[col] = R_all
    
    df_ph_coh_decays[k] = ph_coh_df_loc.mean(axis=1)
    plt.plot(time, ph_coh_df_loc.mean(axis=1), linewidth=5, label=str(k))
plt.xlabel('Time [hours]')
plt.ylabel('Phase coherence of phase differences') #or only phase coherence for circadian
plt.legend(loc='best')
plt.savefig(output+'Phase_coh_ph_diff_decay_entire_image.svg')
plt.show()
        
df_ph_coh_decays.to_excel(path2+'Weighted_ph_coh_ph_dif_entire_image.xlsx')


#%%Phase differences
df_ph_coh_decays = pd.read_excel(path2+'Weighted_ph_coh_ph_dif_entire_image.xlsx', index_col=0, header=0)

slope = []
error = []
rel_final_ph_coh = []
decays = [0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 1]

plt.figure(figsize=(12,10))

for col in decays:
    ph_coh = df_ph_coh_decays[col] 
    time = df_ph_coh_decays.index.values*0.5
    params, covariance = curve_fit(func_fit, time, ph_coh, method='lm',  maxfev=20000)   
    y_fit = func_fit(time, *params)
    plt.plot(time, df_ph_coh_decays[col], linewidth=5, label=str(col))
    plt.plot(time, y_fit, '--', linewidth=5)
    rel_final_ph_coh.append(ph_coh.iloc[-1]/ph_coh[0])
    slope.append(params[0])
    error.append(np.sqrt(np.diag(covariance))[0])
              
plt.xlabel('Time [hours]')    
plt.ylabel('Phase coherence of phase differences')    
plt.legend(loc='best')
plt.savefig(output+'Phase_coh_ph_diff_decay_entire_image_fitting.svg')
plt.show()

plt.figure(figsize=(12,10))
plt.errorbar(decays, slope, yerr=error, markersize=10)
plt.xscale('log', base=10)
# plt.yscale('log', base=10)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Decay rate')
plt.ylabel('Rate phase decoherence')
# plt.savefig(output+'Rate_decoherence_ph_diff_fitting.svg')
plt.show()

plt.figure(figsize=(12,10))
plt.plot(decays, rel_final_ph_coh, 'o', markersize=10)
plt.xscale('log', base=10)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Decay rate')
plt.ylabel('Relative phase coherence (final/initial)')
plt.savefig(output+'Relative_decoherence_ph_diff_fitting.svg')
plt.show()





