#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:57:48 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
import seaborn as sns

dt = 0.5 # the sampling interval, 0.5hours
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path = '.../Division_matrix/'
path2 = '.../Raw_data/'
output = '.../Fig4/'
dose = ['untreated', '5uM', '10uM']
density = ['high']
channel1 = 'circadian'

colors = {'untreated':'gold', '5uM':'tab:orange', '10uM':'tab:blue'}

#%% Changing inhibitor
imt_all_dose = {}

for i in density:
    df = pd.DataFrame()

    for j in dose:
        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Division matrix
        data = pd.read_excel(path+'divisions_'+str(condition)+r'.xlsx')
        time = data.index.values
        
        #Circadian signal
        data1 = pd.read_csv(path2+condition+channel1+r'_filtered.csv', index_col=0)
        
        imt = {str(num): [] for num in range(2, 6)}
        imt_all = []
        
        for col in data:
            if (data1[col].dropna()).index.values[0] == 0:# and (data1[col].dropna()).index.values[-1]>210:
                division_times = ((np.where(data[col]==1)[0]))
            
                imt_single = []
                if 2 <= len(division_times) <= 6:
                    for ii in range(len(division_times)-1):
                        imt_single.append((division_times[ii+1]-division_times[ii])*0.5)
                        imt_all.append((division_times[ii+1]-division_times[ii])*0.5)
                    imt[str(len(division_times))].append(np.mean(imt_single))
        
        imt_all_dose[j] = imt_all
        
        # df_imt = pd.DataFrame.from_dict(imt, orient='index')
        # df_imt = df_imt.transpose()
        # df_imt.columns = [str(num)+' divisions' for num in range(2, 6)]
    
        # fig = plt.figure(figsize=(10,10))
        # df_imt.boxplot(showfliers=False)
        # plt.ylabel('IMT [hours]')
        # # plt.savefig(output+'IMT_boxplot_numdivisions_high_density.svg')
        # plt.show()
        
    max_imt = np.max([max(imt_all_dose[j]) for j in dose])    
    bin_width = max_imt/40
    
    plt.figure(figsize=(10,8))
    for ii in dose:
        print(len(imt_all_dose[ii]))
        print(np.mean(imt_all_dose[ii]), np.std(imt_all_dose[ii]))
        sns.histplot(imt_all_dose[ii], kde=True, stat='density', label= str(ii), color=colors[ii], bins=np.arange(16, max_imt + bin_width, bin_width))
    plt.xlabel('IMT [hours]')
    plt.legend(loc='best')
    # plt.savefig(output+'IMT_distributions_high_density.svg')
    plt.show()
        

#%% Changing density
            
density = ['high','medium','low']
        
imt_all_density = {}

for i in density:
    df = pd.DataFrame()
    
    condition = str(i)+'_density_'+str(dose[0])+'_'
    print(condition)
    
    #Division matrix
    data = pd.read_excel(path+'divisions_'+str(condition)+r'.xlsx')
    time = data.index.values
    
    #Circadian signal
    data1 = pd.read_csv(path2+condition+channel1+r'_filtered.csv', index_col=0)
    
    imt = {str(num): [] for num in range(2, 6)}
    imt_all = []
    
    for col in data:
        if (data1[col].dropna()).index.values[0] == 0:#and (data1[col].dropna()).index.values[-1]>210:
            division_times = ((np.where(data[col]==1)[0]))
        
            imt_single = []
            if 2 <= len(division_times) <= 6:
                for ii in range(len(division_times)-1):
                    imt_single.append((division_times[ii+1]-division_times[ii])*0.5)
                    imt_all.append((division_times[ii+1]-division_times[ii])*0.5)
                imt[str(len(division_times))].append(np.mean(imt_single))
    
    imt_all_density[i] = imt_all
# print(imt_all_density)  
 
colors = {'high':'gold', 'medium':'tab:orange', 'low':'tab:blue'}
       
max_imt = np.max([max(imt_all_density[i]) for i in density])    
bin_width = max_imt/40

plt.figure(figsize=(10,8))
for ii in density:
    print(len(imt_all_density[ii]))
    print(np.mean(imt_all_density[ii]), np.std(imt_all_density[ii]))
    sns.histplot(imt_all_density[ii], kde=True, stat='density', label= str(ii), color=colors[ii], bins=np.arange(16, max_imt + bin_width, bin_width))
plt.xlabel('IMT [hours]')
plt.legend(loc='best')
# plt.savefig(output+'IMT_distributions_changing_density.svg')
plt.show()
               