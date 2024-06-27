#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:53:06 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from pyboat import WAnalyzer
    
plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'      

def exponential_func(t, N0, r):
    K = max(ydata)
    return K/(1+((K-N0)/N0)*np.exp(-r*t))


dt = 1.5 # the sampling interval, 0.5hours
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

path = 'Raw_data/'
output = '...'
file1 = 'U2OS_3osc_dko_1stplate.xlsx'
file2 = 'U2OS_3osc_dko_3rdplate.xlsx'

data1_confl = pd.read_excel(path+file1, sheet_name=0, skiprows=[0], header=[0])
data1_green = pd.read_excel(path+file1, sheet_name=1, skiprows=[0], header=[0]) 

data2_confl = pd.read_excel(path+file2, sheet_name=0, skiprows=[0], header=[0])
data2_green = pd.read_excel(path+file2, sheet_name=1, skiprows=[0], header=[0]) 

#%%Average confluence results 1st plate
data1_confl = data1_confl.drop(columns=['Date Time', 'Elapsed'])

index = (max(data1_confl.idxmin()))

data1_confl_smoothed = pd.DataFrame()

for col in data1_confl:
    ydata = data1_confl[col][index:]
    data1_confl_smoothed[col] = savgol_filter(ydata/ydata[index], window_length=9, polyorder=2)

time = data1_confl_smoothed.index.values*1.5+index
reset = 57#-index*1.5
    
aver_dKo_high_dmso = data1_confl_smoothed[['A1','A2','A3']].mean(axis=1) #'A2', 'A3'
# aver_dKo_low_dmso = data1_confl_smoothed[['B1', 'B2', 'B3']].mean(axis=1)
aver_dKo_high_inh = data1_confl_smoothed[['A4']].mean(axis=1) #, 'A5', 'A6'
# aver_dKo_low_inh = data1_confl_smoothed[['B4', 'B5']].mean(axis=1) #B6

std_dKo_high_dmso = data1_confl_smoothed[['A1', 'A2','A3']].std(axis=1) #'A2', 'A3'
# std_dKo_low_dmso = data1_confl_smoothed[['B1', 'B2', 'B3']].std(axis=1)
std_dKo_high_inh = data1_confl_smoothed[['A4','A5','A6']].std(axis=1) #, 'A5', 'A6'
# std_dKo_low_inh = data1_confl_smoothed[['B4', 'B5']].std(axis=1) #B6

aver_3osc_high_dmso = data1_confl_smoothed[['C1', 'C2', 'C3']].mean(axis=1)
# aver_3osc_low_dmso = data1_confl_smoothed[['D1', 'D2', 'D3']].mean(axis=1)
aver_3osc_high_inh = data1_confl_smoothed[['C5', 'C6']].mean(axis=1) #C4
# aver_3osc_low_inh = data1_confl_smoothed[['D4', 'D5', 'D6']].mean(axis=1)

std_3osc_high_dmso = data1_confl_smoothed[['C1', 'C2','C3']].std(axis=1) 
# std_3osc_low_dmso = data1_confl_smoothed[['D1', 'D2', 'D3']].std(axis=1)
std_3osc_high_inh = data1_confl_smoothed[['C5','C6']].std(axis=1) 
# std_3osc_low_inh = data1_confl_smoothed[['D4', 'D5']].std(axis=1) 

plt.figure(figsize=(10,8))

plt.plot(time, aver_dKo_high_dmso, linewidth=5, color='blue', label='dKO high DMSO')
plt.errorbar(time, aver_dKo_high_dmso, yerr=std_dKo_high_dmso, color='blue', alpha=0.3)
# ydata = aver_dKo_high_dmso
# params, covariance = curve_fit(exponential_func, time, ydata)
# y_fit = exponential_func(time, *params)
# print('r', params[1])
# errors = np.sqrt(np.diag(covariance))
# print('error', errors[1])
# plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)

plt.plot(time, aver_dKo_high_inh, '--', linewidth=5, color='tab:green', label='dKO high inh')
plt.errorbar(time, aver_dKo_high_inh, yerr=std_dKo_high_inh, color='tab:green', alpha=0.3)
# ydata = aver_dKo_high_inh
# params, covariance = curve_fit(exponential_func, time, aver_dKo_high_inh)
# y_fit = exponential_func(time, *params)
# print('r', params[1])
# errors = np.sqrt(np.diag(covariance))
# print('error', errors[1])
# plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)

plt.plot(time, aver_3osc_high_dmso, linewidth=5, color='brown', label='wt high DMSO')
plt.errorbar(time, aver_3osc_high_dmso, yerr=std_3osc_high_dmso, color='brown', alpha=0.3)
# ydata = aver_3osc_high_dmso
# params, covariance = curve_fit(exponential_func, time, ydata)
# y_fit = exponential_func(time, *params)
# print('r', params[1])
# errors = np.sqrt(np.diag(covariance))
# print('error', errors[1])
# plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)

plt.plot(time, aver_3osc_high_inh, '--', linewidth=5, color='tab:orange', label='wt high inh')
plt.errorbar(time, aver_3osc_high_inh, yerr=std_3osc_high_inh, color='tab:orange', alpha=0.3)
# ydata = aver_3osc_high_inh
# params, covariance = curve_fit(exponential_func, time, ydata)
# y_fit = exponential_func(time, *params)
# print('r', params[1])
# errors = np.sqrt(np.diag(covariance))
# print('error', errors[1])
# plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)

plt.xlabel('Time [hours]')
plt.ylabel('Normalized confluence')
plt.legend(loc='best')
plt.xticks([24,48,72,96,120,148])
plt.savefig(output+'Raw_confl_3osc_dko_DMSO_inh_high_dens.svg')
plt.show()


plt.figure(figsize=(10,8))

plt.plot(time, aver_dKo_high_dmso/aver_dKo_high_inh, linewidth=5, label='dKO high')
ydata = aver_dKo_high_dmso/aver_dKo_high_inh
params, covariance = curve_fit(exponential_func, time, ydata)
y_fit = exponential_func(time, *params)
plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)
errors = np.sqrt(np.diag(covariance))
print(errors[1])

plt.plot(time, aver_3osc_high_dmso/aver_3osc_high_inh, linewidth=5, label='3osc high')
ydata = aver_3osc_high_dmso/aver_3osc_high_inh
params, covariance = curve_fit(exponential_func, time, ydata)
y_fit = exponential_func(time, *params)
plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)
errors = np.sqrt(np.diag(covariance))
print(errors[1])

plt.xlabel('Time [hours]')
plt.ylabel('Relative confluence')
plt.legend(loc='best')
plt.savefig(output+'Relative_confl_3osc_dko_DMSO_inh_high_dens_reset0.svg')
plt.show()