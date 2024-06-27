#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:17:38 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em

plt.rcParams.update({'font.size': 26})
plt.rcParams['svg.fonttype'] = 'none'

path = '...Raw_data/'
output = '.../'
dose = ['0uM']
density = ['high']
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
        phases_circ = pd.DataFrame(index=data1.index.values)
        phases_diff = pd.DataFrame(index=data1.index.values)
        
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
                phases_circ[column] = rd['phase']

                
                wAn.compute_spectrum(signal, do_plot=False)
                rd2 = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)
                rd2.set_index(signal.index, inplace=True)
                ridge_results2[column] = rd2
                rd2['traceId'] = column
                df2 = df2.append(rd2)
                
                diff_rel = (rd['phase'] - rd2['phase'])
                diff = np.arctan2(np.sin(diff_rel), np.cos(diff_rel))

                phases_diff[column] = diff
        
        powers_series = em.average_power_distribution(ridge_results1.values(), signal_ids = ridge_results1.keys())
        powers_series = powers_series.reset_index()
    
        
        #%% Periods ---------------------------
        df_period_diff = pd.DataFrame(index=data1.index.values)
        for col in data1:
            period_circ = ridge_results1[col]['periods']
            period_cc = ridge_results2[col]['periods']
            df_period_diff[col] = period_circ-period_cc


        #%% Phase coherence -----------------------------
        phases_circ = phases_circ.T
        phases_diff = phases_diff.T
        R_circ = np.zeros(len(time))
        R_diff = np.zeros(len(time))
                
        for col in phases_circ:
            sum_phi_circ = 0
            sum_phi_all = 0
            count = 0
            for ii in range(len(phases_circ[col])):
                if str(phases_circ[col][ii]) != 'nan':
                    count+=1
                    sum_phi_circ = sum_phi_circ+np.exp(phases_circ[col][ii]*1j)
                    sum_phi_all = sum_phi_all+np.exp(phases_diff[col][ii]*1j)
            R_circ[col] = np.abs(sum_phi_circ/count)
            R_diff[col] = np.abs(sum_phi_all/count)
        R_circ[len(R_circ)-1] = R_circ[len(R_circ)-2]
        R_diff[len(R_diff)-1] = R_diff[len(R_diff)-2]
        

        #%%Overlapped plot --------------------------------

        max_value = max(R_circ)

        tp = 0.5*int(len(time)/2)#0.5*min(range(len(list(R_circ))), key=lambda i: abs(list(R_circ)[i]-max_value/2))
        print(tp)
        
        detuning = wAn.sinc_smooth(df_period_diff.median(axis=1), T_c=20)
        
        fig, ax1 = plt.subplots(figsize=(14, 12))
        ax1.plot(time, R_circ/max(R_circ), linewidth=5, color='tab:brown', label='Circadian clocks')
        ax1.plot(time, R_diff/max(R_diff), linewidth=5, color='tab:purple', label='Phase differences')
        ax1.axvline(x=tp, color='grey', linestyle='--')
        ax1.set_ylabel('Relative phase coherence')
        ax2 = ax1.twinx()
        ax2.plot(time, detuning/(-min(detuning)), linewidth=5, color='tab:grey', label='Detuning circadian cell cycle')
        ax2.set_ylabel('Relative detuning [hours]')
        ax1.set_xlabel('Time [hours]')
        ax1.set_xticks([0,24,48,72,96,120])
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='upper left')
        plt.savefig(output+'Median_period_ph_coh_circ_ph_diff_weak_coupling_experimental.svg')
        plt.show()   
        
        
        #%%Split the experiment in two parts and choose the first/second part of the data
        
        #1st part
        select_index = powers_series.loc[powers_series[0]>10]['index']
        circadian = df1.loc[df1['traceId'].isin(select_index)]
        cell_cycle = df2.loc[df2['traceId'].isin(select_index)]
        
        select_index_2=cell_cycle.loc[cell_cycle['time']<tp]['time'] # change ><60
        cell_cycle=cell_cycle.loc[cell_cycle['time'].isin(select_index_2)]
        circadian=circadian.loc[circadian['time'].isin(select_index_2)]
                
        ph_cell = (cell_cycle['phase']/(2*np.pi)).to_numpy()
        ph_circ = (circadian['phase']/(2*np.pi)).to_numpy()

        fig = plt.figure(figsize = (12,10))
        hist2D = plt.hist2d(ph_circ, ph_cell, bins=np.arange(0,1.01,0.05), cmap='plasma', density=True, vmin=0.5, vmax=1.5)
        plt.colorbar()
        plt.xlabel(r'Circadian phase ($\theta$/2$\pi$)')
        plt.ylabel(r'Cell cycle phase ($\theta$/2$\pi$)')
        # plt.savefig(output+'Phase_lock'+str(condition)+'1st.svg')
        plt.show()
          
        #2nd part             
        select_index = powers_series.loc[powers_series[0]>10]['index']
        circadian = df1.loc[df1['traceId'].isin(select_index)]
        cell_cycle = df2.loc[df2['traceId'].isin(select_index)]

        select_index_2=cell_cycle.loc[cell_cycle['time']>tp]['time'] # change ><60
        cell_cycle=cell_cycle.loc[cell_cycle['time'].isin(select_index_2)]
        circadian=circadian.loc[circadian['time'].isin(select_index_2)]
                
        ph_cell = (cell_cycle['phase']/(2*np.pi)).to_numpy()
        ph_circ = (circadian['phase']/(2*np.pi)).to_numpy()

        fig = plt.figure(figsize = (12,10))
        hist2D = plt.hist2d(ph_circ, ph_cell, bins=np.arange(0,1.01,0.05), cmap='plasma', density=True, vmin=0.5, vmax=1.5)
        plt.colorbar()
        plt.xlabel(r'Circadian phase ($\theta$/2$\pi$)')
        plt.ylabel(r'Cell cycle phase ($\theta$/2$\pi$)')
        # plt.savefig(output+'Phase_lock'+str(condition)+'2nd.svg')
        plt.show()
        