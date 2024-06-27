#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:59:26 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
from pyboat import WAnalyzer
import seaborn as sns
from scipy.optimize import curve_fit

def exponential_func(t, N0, r):
    K = max(df[j])
    return K/(1+((K-N0)/N0)*np.exp(-r*t))

def linear_func(t, a, b):
    return a*t+b

dt = 0.5 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path = '.../Division_matrix/'
output = '.../Fig4/'
dose = ['untreated','5uM','10uM']
density = ['high']
channel2 = 'cell_cycle'

#%%One density
for i in density:
    
    df = pd.DataFrame()
    fit_params = []
    cov_error = []
    
    fit_params_l = []
    cov_error_l = []

    plt.figure(figsize=(12,10))

    for j in dose:
        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Division matrix
        data = pd.read_excel(path+'divisions_'+str(condition)+r'.xlsx')
        time = data.index.values
        print(len(data.columns))
            
        median_distr = pd.DataFrame()
        for a in range(100):
            num_Ids = random.sample(list(data.columns), 750) #high=750, medium=112, low=49
            time_points = len(time)
            division_profile = np.zeros((len(num_Ids),time_points))
            
            division_profile = np.array([data[col].values for col in num_Ids])
                            
            distr = (np.sum(division_profile, axis=0))
            median_distr[a] = np.cumsum(distr)

        average_distr = median_distr.mean(axis=1)
    
        df[j] = average_distr   
        maximum = (df.values.max())
        
        params_l, covariance_l = curve_fit(linear_func, time[0:150], df[j][0:150])
        errors_l = np.sqrt(np.diag(covariance_l))
        y_fit_l = linear_func(time[0:150], *params_l)
        fit_params_l.append(params_l[0])
        cov_error_l.append(errors_l[0])
    
        params, covariance = curve_fit(exponential_func, time, df[j])
        errors = np.sqrt(np.diag(covariance))
        y_fit = exponential_func(time, *params)
        fit_params.append(params[1])
        cov_error.append(errors[1])
        print(errors[1])

        plt.plot(time*0.5, df[j]/maximum, label = str(j), linewidth=4)
        plt.plot(time[24:len(time)]*0.5, y_fit[24:len(time)]/maximum, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)
        # plt.plot(time[0:150]*0.5, y_fit_l, '--', label='linear fit: {:.3f}'.format(params_l[0]), linewidth=5, alpha=0.5)
    plt.xticks([0,24,48,72,96,120])
    plt.xlabel('Time(h)')
    plt.ylabel('Cumulative distribution of division events')
    plt.legend(loc = 'best')
    # plt.savefig(output+'Cumul_distr_all_divisions_0_5_10uM.svg')
    plt.show()
    
    print(fit_params)
    print(cov_error)
    
    plt.figure(figsize=(10,8))
    plt.errorbar(dose, fit_params, yerr=cov_error, fmt='o')
    plt.xlabel('Inhibitor')
    plt.ylabel('Rate of the logistic function')
    # plt.savefig(output+'Logistic_fit_growth_curves.svg')
    plt.show()
        
    # plt.figure(figsize=(10,8))
    # plt.errorbar(dose, fit_params_l, yerr=cov_error_l, fmt='o')
    # plt.xlabel('Inhibitor')
    # plt.ylabel('Growth rate')
    # # plt.savefig(output+'Growth_rates_curves.svg')
    # plt.show()
    
    plt.figure(figsize=(12,10))    
    for col in df:    
        plt.plot(time*0.5, wAn.sinc_smooth(wAn.sinc_detrend(df[col]/maximum, T_c=40), T_c=10), label = str(col), linewidth=4)
    plt.xlabel('Time(h)')
    plt.ylabel('Detrended cumulative distribution of division events')
    plt.xticks([0,24,48,72,96,120])
    plt.legend(loc = 'best')
    # plt.savefig(output+'Detrended_Cumul_distr_all_divisions_0_5_10uM.svg')
    plt.show()
    
    periods = {}
    
    for col in df:    
        signal = wAn.sinc_detrend(df[col]/maximum, T_c=42)
        wAn.compute_spectrum(signal, do_plot=False, draw_coi=False)
        rd = wAn.get_maxRidge(power_thresh = 0, smoothing_wsize=20)
        # wAn.draw_Ridge()
        periods[col] = rd['periods']
        
    plt.figure(figsize=(12,10)) 
    for col in df:
        plt.plot(periods[col].index.values*0.5, periods[col]/periods['untreated'], label = str(col), linewidth=4)
    plt.xlim([-5,65])
    plt.ylabel('Period detrended cumulative division events [hours]')
    plt.legend(loc = 'best')
    # plt.savefig(output+'Period_detrended_Cumul_distr_all_divisions_0_5_10uM.svg')
    plt.show()
        
    periods_df = pd.DataFrame().from_dict(periods, orient='index')
    periods_df = periods_df.T
    
    for col in periods_df:
        print(col, np.mean(periods_df[col]), np.std(periods_df[col]))

    plt.figure(figsize=(12,10)) 
    sns.boxplot(data=periods_df, showfliers=False)
    # for col in df:
    #     plt.plot(periods[col].index.values*0.5, periods[col], label = str(col), linewidth=4)
    # plt.ylim([22,35])
    plt.ylabel('Period detrended cumulative division events [hours]')
    # plt.legend(loc = 'best')
    # plt.savefig(output+'Boxplot_period_detrended_cumul_distr_cut_yaxos.svg')
    plt.show()
    
    
#%%Density

density = ['high','medium','low']

plt.figure(figsize=(12,10))

df = pd.DataFrame()
fit_params = []
cov_error = []

fit_params_l = []
cov_error_l = []

for j in density:
    condition = str(j)+'_density_untreated_'
    print(condition)
    
    #Division matrix
    data = pd.read_excel(path+'divisions_'+str(condition)+r'.xlsx')
    time = data.index.values
    print(len(data.columns))
        
    median_distr = pd.DataFrame()
    for a in range(100):
        num_Ids = random.sample(list(data.columns), 89) #high=750, medium=112, low=49
        time_points = len(time)
        division_profile = np.zeros((len(num_Ids),time_points))
        
        division_profile = np.array([data[col].values for col in num_Ids])
                        
        distr = (np.sum(division_profile, axis=0))
        median_distr[a] = np.cumsum(distr)

    average_distr = median_distr.mean(axis=1)

    df[j] = average_distr   
    maximum = (df.values.max())
    
    params_l, covariance_l = curve_fit(linear_func, time[0:150], df[j][0:150])
    errors_l = np.sqrt(np.diag(covariance_l))
    y_fit_l = linear_func(time[0:150], *params_l)
    fit_params_l.append(params_l[0])
    cov_error_l.append(errors_l[0])

    params, covariance = curve_fit(exponential_func, time, df[j])
    errors = np.sqrt(np.diag(covariance))
    y_fit = exponential_func(time, *params)
    fit_params.append(params[1])
    cov_error.append(errors[1])
    print(errors[1])


    plt.plot(time*0.5, df[j]/maximum, label = str(j), linewidth=4)
    plt.plot(time[24:len(time)]*0.5, y_fit[24:len(time)]/maximum, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, alpha=0.5)
    # plt.plot(time[0:150]*0.5, y_fit_l, '--', label='linear fit: {:.3f}'.format(params_l[0]), linewidth=5, alpha=0.5)

print(fit_params)
print(cov_error)

plt.xticks([0,24,48,72,96,120])
plt.xlabel('Time(h)')
plt.ylabel('Cumulative distribution of division events')
plt.legend(loc = 'best')
# plt.savefig(output+'Cumul_distr_all_divisions_change_density.svg')
plt.show()

plt.figure(figsize=(10,8))
plt.errorbar(dose, fit_params, yerr=cov_error, fmt='o')
plt.xlabel('Inhibitor')
plt.ylabel('Rate of the logistic function')
# plt.savefig(output+'Logistic_fit_growth_curves_density.svg')
plt.show()
    
# plt.figure(figsize=(10,8))
# plt.errorbar(dose, fit_params_l, yerr=cov_error_l, fmt='o')
# plt.xlabel('Inhibitor')
# plt.ylabel('Growth rate')
# # plt.savefig(output+'Growth_rates_curves.svg')
# plt.show()

plt.figure(figsize=(12,10))    
for col in df:    
    plt.plot(time*0.5, wAn.sinc_smooth(wAn.sinc_detrend(df[col]/maximum, T_c=40), T_c=10), label = str(col), linewidth=4)
plt.xlabel('Time(h)')
plt.ylabel('Detrended cumulative distribution of division events')
plt.xticks([0,24,48,72,96,120])
plt.legend(loc = 'best')
# plt.savefig(output+'Detrended_Cumul_distr_change_density.svg')
plt.show()

periods = {}

for col in df:    
    signal = wAn.sinc_detrend(df[col]/maximum, T_c=42)
    wAn.compute_spectrum(signal, do_plot=False)
    rd = wAn.get_maxRidge(power_thresh = 0, smoothing_wsize=20)
    # wAn.draw_Ridge()
    periods[col] = rd['periods']
    
# plt.figure(figsize=(12,10)) 
# for col in df:
#     plt.plot(periods[col].index.values*0.5, periods[col]/periods['untreated'], label = str(col), linewidth=4)
# plt.xlim([-5,65])
# plt.ylabel('Period detrended cumulative division events [hours]')
# plt.legend(loc = 'best')
# # plt.savefig(output+'Period_detrended_Cumul_distr_all_divisions_0_5_10uM.svg')
# plt.show()
    
periods_df = pd.DataFrame().from_dict(periods, orient='index')
periods_df = periods_df.T

for col in periods_df:
    print(col, np.mean(periods_df[col]), np.std(periods_df[col]))

plt.figure(figsize=(12,10)) 
sns.boxplot(data=periods_df, showfliers=False)
# for col in df:
#     plt.plot(periods[col].index.values*0.5, periods[col], label = str(col), linewidth=4)
# plt.ylim([22,35])
plt.ylabel('Period detrended cumulative division events [hours]')
# plt.legend(loc = 'best')
# plt.savefig(output+'Boxplot_period_detrended_cumul_distr_change_density.svg')
plt.show()