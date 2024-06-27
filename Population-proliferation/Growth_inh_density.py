#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:46:51 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from pyboat import WAnalyzer
import seaborn as sns
    
plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'    
  
def exponential_func(t, N0, r):
    K = max(ydata_smoothed)
    return K/(1+((K-N0)/N0)*np.exp(-r*t))

def linear_func(t, a, b):
    return a*t+b

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def exp_func(b, x):
    return 0.068*np.exp(b*x)

def exp_func2(b, x):
    return 0.13*np.exp(-b*x)

dt = 1 
lowT = 16
highT = 40
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')
            
#%%Data
path = 'Raw_data/'
output = '.../'

file = 'Circadian_extra_coupling.xlsx'

data1_confl = pd.read_excel(path+file, sheet_name=0, skiprows=[0], header=[0])
data1_green = pd.read_excel(path+file, sheet_name=1, skiprows=[0], header=[0]) 

data2_confl = pd.read_excel(path+file, sheet_name=4, skiprows=[0], header=[0])
data2_green = pd.read_excel(path+file, sheet_name=5, skiprows=[0], header=[0]) 

#%%Average results for confluence
data1_confl = data1_confl.drop(columns=['Date Time', 'Elapsed'])
data2_confl = data2_confl.drop(columns=['Date Time', 'Elapsed'])
data_confl = (data1_confl+data2_confl)/2
data_confl = data_confl.drop([2])

line_colors = ['tab:green', 'tab:orange', 'tab:blue']
labels1 = ['10uM', '5uM', 'untreated']
labels2 = ['high', 'medium', 'low']

columns1 = ['A1', 'A2', 'A6']
columns2 = ['A6', 'B6', 'C6']

fit_params = []
cov_error = []

fit_params_l = []
cov_error_l= []

doubling_times = []

plt.figure(figsize=(10,8))
for col, label, color in zip(columns2, labels2, line_colors):
    minimum = data_confl[col].min()
    index = data_confl[col].idxmin()
    ydata = data_confl[col][index:]/minimum
    time = ydata.index.values

    ydata_smoothed = savgol_filter(ydata, window_length=9, polyorder=2)
    params, covariance = curve_fit(exponential_func, time, ydata_smoothed)
    print(params)
    errors = np.sqrt(np.diag(covariance))
    y_fit = exponential_func(time, *params)
    fit_params.append(params[1])
    cov_error.append(errors[1])
    doubling_times.append(np.log(2)/params[1])
    
    plt.plot(time, ydata_smoothed, label=label, linewidth=5, color=color)
    plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, color=color, alpha=0.5)

plt.ylabel('Confluence (%)')
plt.xlabel('Time [hours]')
plt.xticks([0,24,48,72,96,120])
plt.legend(loc='best')
# plt.savefig(output+'Confluence_change_density.svg')
plt.show()

print(fit_params)
print(cov_error)

plt.figure(figsize=(10,8))
plt.errorbar(labels2, fit_params/fit_params[-1], yerr=cov_error, fmt='o')
# params, covariance = curve_fit(exp_func2, [0,1,2], fit_params) ##change exp_func or exp_func2
# print(params)
# error = np.sqrt(np.diag(covariance))
# y_fit = exp_func2(*params, np.linspace(0, 2, 100))
# plt.plot(np.linspace(0, 2, 100), y_fit, '--', label=r'$e^{{(-{:.2f} \pm {:.2f}) \cdot r}}$'.format(params[0], error[0]))
plt.xlabel('Seeding density') #Inhibitor concentration
plt.ylabel('Growth rate')
# plt.yscale('log')
plt.legend(loc='best')
# plt.savefig(output+'Logistic_fit_growth_curves_density.svg')
plt.show()

# plt.figure(figsize=(10,8))
# plt.plot(labels2, doubling_times, 'o', markersize=15)
# plt.xlabel('Seeding density')
# plt.ylabel('Doubling times [hours]')
# # plt.savefig(output+'Doubling_times_density.svg')
# plt.show()

# plt.figure(figsize=(10,8))
# plt.errorbar(labels1, fit_params_l, yerr=cov_error_l, fmt='o')
# plt.xlabel('Inhibitor')
# plt.ylabel('Linear slope')
# # plt.savefig(output+'Growth_rates_curves_incucyte_exp.svg')
# plt.show()


#%%Average results for circadian signal
data1_green = data1_green.drop(columns=['Date Time', 'Elapsed'])
data2_green = data2_green.drop(columns=['Date Time', 'Elapsed'])
data_green = (data1_green+data2_green)/2

# line_colors = ['tab:green', 'tab:orange', 'tab:blue']
labels1 = ['10uM', '5uM', 'untreated']
labels2 = ['high', 'medium', 'low']
line_colors = ['tab:green', 'tab:orange', 'tab:blue']

columns1 = ['A1', 'A2', 'A6']
columns2 = ['A6', 'B6', 'C6']

plt.figure(figsize=(10, 8))
for col, label, color in zip(columns1, labels1, line_colors):
    minimum = data_confl[col].min()
    index = data_green[col].idxmin()
    ydata = data_green[col][index:]
    time = ydata.index.values

    ydata_smoothed = savgol_filter(ydata, window_length=9, polyorder=2)

    plt.plot(time, ydata_smoothed/ydata_smoothed[0], label=label, linewidth=5, color=color)

plt.ylabel('Circadian signal [a.u.]')
plt.xlabel('Time [hours]')
plt.xticks([0,24,48,72,96,120])
plt.legend(loc='best')
# plt.savefig(output+'Circadian_inhibitor.svg')
plt.show()

periods_df = pd.DataFrame()
amplitudes_df = pd.DataFrame()
for col in (columns1):
    ydata = data_green[col]
    ydata_smoothed = savgol_filter(ydata, window_length=9, polyorder=2)
    ydata_det = wAn.sinc_detrend(ydata_smoothed, T_c=50)
    ydata_norm = wAn.normalize_amplitude(ydata_det, window_size=70)
    
    wAn.compute_spectrum(ydata_norm, do_plot=False)
    wAn.get_maxRidge(power_thresh = 5, smoothing_wsize=15)
    # wAn.draw_Ridge()
    # plt.savefig(output+'Power_spectrum_'+str(col)+'.svg')
    
    rd = wAn.ridge_data
    periods_df[col] = rd['periods']
    amplitudes_df[col] = rd['amplitude']

plt.figure(figsize=(10, 8))
for col, label, color in zip(columns1, labels1, line_colors):
    time = periods_df[col].index.values
    plt.plot(time, periods_df[col], label=label, linewidth=5, color=color)
plt.ylabel('Period [hours]')
plt.xlabel('Time [hours]')
plt.legend(loc='best')
# plt.savefig(output+'Periods_inhibitor.svg')
plt.show()

    
plt.figure(figsize=(10,8))
for col, label in zip(columns2, labels2):
    plt.plot(periods_df[col], amplitudes_df[col], 'o', label=label, alpha=0.5) # c=vect, cmap=colormaps[kk],   
plt.xlabel('Periods [hours]')
plt.ylabel('Amplitude')
plt.legend(loc='best')
# plt.savefig(output+'Period_amplitude_changing_density.svg')
plt.show()

periods_df.columns = labels1
plt.figure(figsize=(10, 8))
sns.boxplot(data=periods_df)
plt.ylabel('Period [hours]')
plt.xlabel('Density')
# plt.savefig(output+'Periods_changing_inhibitor_boxplot.svg')
plt.show()