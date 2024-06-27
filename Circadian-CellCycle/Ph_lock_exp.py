#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:58:03 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em

plt.rcParams.update({'font.size': 20})
plt.rcParams['svg.fonttype'] = 'none'

path = '.../Raw_data/'
output = '.../Fig2/'
dose = ['untreated','5uM','10uM'] 
density = ['high','medium','low']
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
        time = data1.index.values*0.5
        
        ridge_results1 = {}
        df1 = pd.DataFrame()
        ridge_results2 = {}
        df2 = pd.DataFrame()

        for column in data1:
            norm_signal = data1[column].dropna()
            signal = data2[column].dropna()

            if len(norm_signal) > 100:
                wAn.compute_spectrum(norm_signal, do_plot=False)
                rd = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)
                rd.set_index(norm_signal.index, inplace=True)
                ridge_results1[column] = rd
                rd['traceId'] = column
                df1 = df1.append(rd)
                
                wAn.compute_spectrum(signal, do_plot=False)
                rd2 = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)
                rd2.set_index(signal.index, inplace=True)
                ridge_results2[column] = rd2
                rd2['traceId'] = column
                df2 = df2.append(rd2)
        
        powers_series = em.average_power_distribution(ridge_results1.values(), signal_ids = ridge_results1.keys())
        powers_series = powers_series.reset_index()
        
        select_index = powers_series.loc[powers_series[0]>10]['index']
        circadian = df1.loc[df1['traceId'].isin(select_index)]
        cell_cycle = df2.loc[df2['traceId'].isin(select_index)]

        # #Split the experiment in two parts and choose the first/second part of the data
        # select_index_2=cell_cycle.loc[cell_cycle['time']>60]['time'] # change ><60
        # cell_cycle=cell_cycle.loc[cell_cycle['time'].isin(select_index_2)]
        # circadian=circadian.loc[circadian['time'].isin(select_index_2)]
                
        #Normalization
        ph_cell = (cell_cycle['phase']/(2*np.pi)).to_numpy()
        ph_circ = (circadian['phase']/(2*np.pi)).to_numpy()

        fig = plt.figure(figsize = (12,10))
        hist2D = plt.hist2d(ph_circ, ph_cell, bins=np.arange(0, 1.01, 0.05), cmap='plasma', density=True, vmin=0.5, vmax=1.5)
        plt.colorbar()
        plt.xlabel(r'Circadian phase ($\theta$/2$\pi$)')
        plt.ylabel(r'Cell cycle phase ($\theta$/2$\pi$)')
        plt.savefig(output+'Phase_lock'+str(condition)+'.svg')
        plt.show()
        
