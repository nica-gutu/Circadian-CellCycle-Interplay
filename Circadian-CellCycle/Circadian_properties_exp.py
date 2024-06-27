#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:52:13 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em
import seaborn as sns

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path = '...Raw_data/'
output = '.../Fig2/'
dose = ['untreated', '5uM', '10uM']
density = ['high']
channel1 = 'circadian'
channel2 = 'cell_cycle'
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:grey']

dt = 0.5 
lowT = 10
highT = 50
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

a = 3 #column of ridge results
magn = 'Phase coherence'
value = 'R'

#%%Changing dose within same density

periods_cond = {}
for i in density:
    count=0
    
    fig = plt.figure(figsize = (12,10))
    
    for j in dose:

        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Circadian data
        data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
        #Cell-cycle data
        data2 = pd.read_csv(path+condition+channel2+r'_filtered.csv', index_col=0)
            
        time=data1.index.values*0.5
        
        ridge_results1 = {}
        df1 = pd.DataFrame()
        ridge_results2 = {}
        df2 = pd.DataFrame()
        
        periods_mean = []

        for column in data1:
            norm_signal = data1[column].dropna()
            signal = data2[column].dropna()

            if len(norm_signal)>100:
                wAn.compute_spectrum(norm_signal,do_plot=False)
                rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
                rd.set_index(norm_signal.index,inplace=True)
                ridge_results1[column] = rd
                rd['traceId'] = column
                df1 = df1.append(rd)
                
                for ii in rd['periods'].index.values:
                    if 16 < rd['periods'][ii] < 50:
                        periods_mean.append(rd['periods'][ii])
                
                wAn.compute_spectrum(signal,do_plot=False)
                rd2 = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
                rd2.set_index(signal.index,inplace=True)
                ridge_results2[column] = rd2
                rd2['traceId'] = column
                df2 = data2.append(rd2)
        
        periods_cond[j] = periods_mean
           
        powers_series = em.average_power_distribution(ridge_results1.values(), signal_ids = ridge_results1.keys())
        high_power_ids = powers_series[powers_series > 10].index
        high_power_ridge_results = [ridge_results1[i] for i in high_power_ids]
        res = em.get_ensemble_dynamics(high_power_ridge_results)
        # pl.ensemble_dynamics(*res)
        # plt.show()
        
        print('Sample size: ', len(ridge_results1))

        print('Maximum ph coh: ', np.where(res[a][str(value)]==max(res[a][str(value)]))[0])
        print('Minimmum ph coh: ', np.where(res[a][str(value)]==min(res[a][str(value)]))[0])
        
        plt.plot(time, res[a][str(value)], color=colors[count], label=str(j), linewidth=5)
        # plt.fill_between(time, res[a]['Q1'], res[a]['Q3'], color=colors[count], alpha=.1)
        count+=1
    plt.xlabel('Time [hours]')
    plt.ylabel(str(magn))
    plt.legend(loc='best')
    # plt.savefig(output+'ED_'+str(a)+'_'+str(i)+'0_10uM.svg')
    plt.show()
                       
periods_cond_df = pd.DataFrame().from_dict(periods_cond, orient='index')
periods_cond_df = periods_cond_df.T

print(periods_cond_df)

plt.figure(figsize=(12,10))
sns.boxplot(data=periods_cond_df, showfliers = False)
plt.ylabel('Periods [hours]')
plt.xlabel('Inhibitor dose')
# plt.savefig(output+'Boxplot_instantaneous_periods_high_density.svg')
plt.show()


#%% Changing density        
density = ['high', 'medium', 'low']

fig = plt.figure(figsize = (12,10))
count = 0

periods_cond = {}

for jj in density:
    condition = str(jj)+'_density_'+str(dose[0])+'_'
    print(condition)
    
    #Circadian data
    data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
    #Cell-cycle data
    data2 = pd.read_csv(path+condition+channel2+r'_filtered.csv', index_col=0)
        
    time=data1.index.values*0.5
    
    ridge_results1 = {}
    df1 = pd.DataFrame()
    ridge_results2 = {}
    df2 = pd.DataFrame()
        
    periods_mean = []

    for column in data1:
        norm_signal = data1[column].dropna()
        signal = data2[column].dropna()

        if len(norm_signal)>100:
            wAn.compute_spectrum(norm_signal,do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            rd.set_index(norm_signal.index,inplace=True)
            ridge_results1[column] = rd
            rd['traceId'] = column
            df1 = df1.append(rd)
            
            for ii in rd['periods'].index.values:
                if 16 < rd['periods'][ii] < 50:
                    periods_mean.append(rd['periods'][ii])

            wAn.compute_spectrum(signal,do_plot=False)
            rd2 = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            rd2.set_index(signal.index,inplace=True)
            ridge_results2[column] = rd2
            rd2['traceId'] = column
            df2 = data2.append(rd2)
    
    periods_cond[jj] = periods_mean
   
    powers_series = em.average_power_distribution(ridge_results1.values(), signal_ids = ridge_results1.keys())
    high_power_ids = powers_series[powers_series > 10].index
    high_power_ridge_results = [ridge_results1[i] for i in high_power_ids]
    res = em.get_ensemble_dynamics(high_power_ridge_results)
    # pl.ensemble_dynamics(*res)
    # plt.show()
    
    print('Sample size: ', len(ridge_results1))
    
    print('Maximum ph coh: ', np.where(res[a][str(value)]==max(res[a][str(value)]))[0])
    print('Minimmum ph coh: ', np.where(res[a][str(value)]==min(res[a][str(value)]))[0])

    plt.plot(time, res[a][str(value)], color=colors[count], label=str(jj), linewidth=5)
    # plt.fill_between(time, res[a]['Q1'], res[a]['Q3'], color=colors[count], alpha=.1)
    count+=1
    
plt.xlabel('Time [hours]')
plt.ylabel(str(magn))
plt.legend(loc='best')
# plt.savefig(output+'ED_'+str(a)+'_changing_density.svg')
plt.show()

periods_cond_df = pd.DataFrame().from_dict(periods_cond, orient='index')
periods_cond_df = periods_cond_df.T

print(periods_cond_df)

colors = {'high': 'tab:orange', 'medium': 'tab:green', 'low': 'tab:blue'}

plt.figure(figsize=(12,10))
sns.boxplot(data=periods_cond_df, showfliers = False)
plt.ylabel('Periods [hours]')
plt.xlabel('Initial density')
# plt.savefig(output+'Boxplot_instantaneous_periods_changing_density.svg')
plt.show()
