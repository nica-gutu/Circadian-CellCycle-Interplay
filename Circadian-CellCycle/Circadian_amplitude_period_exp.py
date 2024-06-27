#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:54:22 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path = '...Raw_data/'
output = '...Fig2/'
dose = ['untreated', '5uM', '10uM']
density = ['high']
channel1 = 'circadian'
channel2 = 'cell_cycle'
colors = ['tab:blue', 'tab:orange', 'tab:green']

dt = 0.5 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#%%Changing dose
for i in density:
    plt.figure(figsize = (12,10))
    periods_all = pd.DataFrame()
    periods_rd = pd.DataFrame()
    count = 0
    for j in dose:

        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Circadian data
        data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
        #Cell-cycle data
        data2 = pd.read_csv(path+condition+channel2+r'_filtered.csv', index_col=0)
            
        time = data1.index.values*0.5
        
        ridge_results1 = {}
        ridge_results2 = {}

        for column in data1:
            norm_signal = data1[column].dropna()
            signal = data2[column].dropna()

            if len(norm_signal)>100:
                wAn.compute_spectrum(norm_signal,do_plot=False)
                rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
                rd.set_index(norm_signal.index,inplace=True)
                ridge_results1[column] = rd
                
                wAn.compute_spectrum(signal,do_plot=False)
                rd2 = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
                rd2.set_index(signal.index,inplace=True)
                ridge_results2[column] = rd2
                
        powers_series = em.average_power_distribution(ridge_results1.values(), signal_ids = ridge_results1.keys())
        high_power_ids = powers_series[powers_series > 10].index        
        
        osc = 0
        for col in high_power_ids:
            osc += data1[col].fillna(0)
            
        signal = wAn.sinc_smooth(osc/len(high_power_ids), T_c=10)
        plt.plot(time, signal, linewidth=5, label=condition)

        wAn.compute_spectrum(signal, do_plot=False)
        rd = wAn.get_maxRidge(power_thresh = 10, smoothing_wsize=50)
        rd_all = wAn.get_maxRidge(power_thresh = 0, smoothing_wsize=50)
        periods_rd[count] = rd['periods']   
        periods_all[count] = rd_all['periods']   
        count+=1

    plt.xlabel('Time [hours]')
    plt.xticks([0,24,48,72,96,120])
    plt.ylabel('Mean circadian signal [a.u.]')
    plt.legend(loc='best')
    plt.savefig(output+'Mean_circ_osc_'+str(i)+'.svg')
    plt.show()
            
    plt.figure(figsize=(10,8))
    for j in range(len(dose)):
        plt.plot(periods_rd[j].index.values*0.5, periods_rd[j], linewidth=5, color=colors[j], label=str(dose[j]))    
        # plt.plot(time, periods_all[j], '--', linewidth=5, color=colors[j], alpha=0.5, label=str(dose[j]))    
    plt.xlabel('Time [hours]')
    plt.xticks([0,24,48,72])#,96,120])
    plt.ylabel('Period')
    plt.legend(loc='best')
    plt.savefig(output+'Period_mean_circ_osc_'+str(i)+'.svg')
    plt.show()

#%% Changing density        
density = ['high', 'medium', 'low']

fig = plt.figure(figsize = (12,10))
count = 0

periods_rd = pd.DataFrame()

for jj in density:
    condition = str(jj)+'_density_untreated_'
    print(condition)
    
    #Circadian data
    data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
        
    time=data1.index.values*0.5
    
    ridge_results1 = {}
    df1 = pd.DataFrame()

    for column in data1:
        norm_signal = data1[column].dropna()

        if len(norm_signal)>100:
            wAn.compute_spectrum(norm_signal,do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            rd.set_index(norm_signal.index,inplace=True)
            ridge_results1[column] = rd
            rd['traceId'] = column
            
    powers_series = em.average_power_distribution(ridge_results1.values(), signal_ids = ridge_results1.keys())
    high_power_ids = powers_series[powers_series > 10].index        
    
    osc = 0
    for col in high_power_ids:
        osc += data1[col].fillna(0)
            
    signal = wAn.sinc_smooth(osc/len(high_power_ids), T_c=10)
    plt.plot(time, signal, linewidth=5, label=str(jj)+' density')

    wAn.compute_spectrum(signal, do_plot=False)
    rd = wAn.get_maxRidge(power_thresh = 10, smoothing_wsize=50)
    rd_all = wAn.get_maxRidge(power_thresh = 0, smoothing_wsize=50)
    periods_rd[count] = rd['periods']   
    count+=1

plt.xlabel('Time [hours]')
plt.xticks([0,24,48,72,96,120])
plt.ylabel('Mean circadian signal [a.u.]')
plt.legend(loc='best')
plt.savefig(output+'Mean_circ_osc_densities.svg')
plt.show()
            
plt.figure(figsize=(10,8))
for j in range(len(density)):
    plt.plot(periods_rd[j].index.values*0.5, periods_rd[j], linewidth=5, color=colors[j], label=str(density[j])+' density')    
plt.xlabel('Time [hours]')
plt.xticks([0,24,48,72])#,96,120])
plt.ylabel('Period')
plt.legend(loc='best')
plt.savefig(output+'Period_mean_circ_osc_densities.svg')
plt.show()
