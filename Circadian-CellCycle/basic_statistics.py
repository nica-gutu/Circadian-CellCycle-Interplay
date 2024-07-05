#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 16:07:04 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
import seaborn as sns

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def calculate_statistics(data):
    mean_val = np.mean(data)
    median_val = np.median(data)
    std_val = np.std(data)
    coeff_var = std_val / mean_val
    return mean_val, median_val, std_val, coeff_var

path = 'Raw_data/'
output = 'Output/'
dose = ['untreated', '5uM','10uM'] #'untreated',, '1.25uM','2.5uM',
density = ['high'] #,'medium','low'
channel1 = 'circadian'

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
        
        amplitudes = []
        periods_list = []
        powers = []

        for column in data1:
            norm_signal = data1[column].dropna()
                                    
            wAn.compute_spectrum(norm_signal, do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            amplitude = np.mean(rd['amplitude'])
            
            # if amplitude>3:
            #     plt.plot(norm_signal.index.values*0.5, norm_signal)
            #     plt.show()

            period = np.mean(rd['periods'])
            power = np.mean(rd['power'])

            amplitudes.append(amplitude)
            periods_list.append(period)
            powers.append(power)           
            
        amplitude_stats = calculate_statistics(amplitudes)
        period_stats = calculate_statistics(periods_list)
        power_stats = calculate_statistics(powers)
        
        # Plot the distributions
        fig, axes = plt.subplots(1, 3, figsize=(22, 8))
        
        # Plot Amplitude
        sns.histplot(amplitudes, bins=50, kde=True, stat='density', ax=axes[0], color='tab:orange', edgecolor='black')
        axes[0].set_xlabel('Detrended mean amplitude [a.u.]')
        axes[0].set_ylabel('Density')
        axes[0].legend([f'Mean: {amplitude_stats[0]:.2f}\nMedian: {amplitude_stats[1]:.2f}\n'
                        f'SD: {amplitude_stats[2]:.2f}\nCV.: {amplitude_stats[3]:.2f}'])
        
        # Plot Period
        sns.histplot(periods_list, bins=30, kde=True, stat='density', ax=axes[1], color='tab:olive', edgecolor='black')
        axes[1].set_xlabel('Mean period [hours]')
        axes[1].legend([f'Mean: {period_stats[0]:.2f}\nMedian: {period_stats[1]:.2f}\n'
                        f'SD: {period_stats[2]:.2f}\nCV.: {period_stats[3]:.2f}'])
        
        # Plot Power
        sns.histplot(powers, bins=30, kde=True, stat='density', ax=axes[2], color='tab:purple', edgecolor='black')
        axes[2].set_xlabel('Mean power [a.u.]')
        axes[2].legend([f'Mean: {power_stats[0]:.2f}\nMedian: {power_stats[1]:.2f}\n'
                        f'SD: {power_stats[2]:.2f}\nCV.: {power_stats[3]:.2f}'])
        
        plt.tight_layout()
        # plt.savefig(output+'Basic_statistics_'+str(j)+'.svg')
        plt.show()        

