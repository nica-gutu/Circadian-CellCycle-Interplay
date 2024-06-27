#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:36:30 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

output = '.../Fig2/'

dt = 0.1 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#Parameters
Nosc = 200
dt = 0.1
tf = 120

gg = 1
ll = 1
a0 = 1
acc = 1

kappa_vect = [0, 0.0005, 0.001, 0.05]
eps = 0.01

periods = pd.DataFrame()
phase_coh = pd.DataFrame()
mean_osc = pd.DataFrame()

plt.figure(figsize=(12,10))

for kk in range(len(kappa_vect)):
    kappa = kappa_vect[kk]
    
    #Vectors initialization
    t = np.arange(0,tf,step=dt)
    X = np.zeros((Nosc,len(t)))
    XX = np.zeros((Nosc,len(t)))
    Y = np.zeros((Nosc,len(t)))
    YY = np.zeros((Nosc,len(t)))
    
    #Initial position
    X[:,0] = np.random.uniform(-0.3,1,size=(Nosc))
    XX[:,0] = np.random.uniform(-1,1,size=(Nosc))
    Y[:,0] = np.random.uniform(-0.5,1,size=(Nosc))
    YY[:,0] = np.random.uniform(-1,1,size=(Nosc))
    
    #Period distribution
    mu1, sigma1 = 22, np.sqrt(6)
    period_circ = np.random.normal(mu1, sigma1, Nosc)
        
    pert_resset = 0

    # Integration over time and population
    for i in range(len(t) - 1):
        sum_x = 0
        sum_y = 0
        
        if int((tf*0.0)/dt)<i<int((tf*0.05)/dt):
            gaussian = 0.01*np.exp(-(i*dt)**2/(2*1**2))
            pert_resset += gaussian
        else:
            pert_resset = 0
                
        for m in range(Nosc):
            middle = -gg*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0)
            X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/period_circ[m]+sum_x*kappa/(2*Nosc)+pert_resset)
            Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/period_circ[m]+sum_y*kappa/(2*Nosc))
                
            for k in range(Nosc):
                sum_x = sum_x+X[k,i] 
                sum_y = sum_y+Y[k,i] 
    osc = 0
    for ii in range(Nosc):
        osc += np.sin(np.arctan2(Y[ii], X[ii]))
    mean_osc[kk] = osc/Nosc
    
    wAn.compute_spectrum(osc/Nosc, do_plot=False)
    rd = wAn.get_maxRidge(power_thresh = 0, smoothing_wsize=20)
    periods[kk] = rd['periods']   
    
    ph = pd.DataFrame(columns=np.arange(0, Nosc, step=1))
    for m in range(Nosc):
        ph_osc = []
        for n in range(len(t)):
            phases = np.arctan2(np.sin(np.arctan2(Y[m,n], X[m,n])), np.cos(np.arctan2(Y[m,n], X[m,n])))
            ph_osc.append(phases)   
        ph[m] = ph_osc
          
    ph = ph.T 
    R_all = np.zeros(len(t))
    for col in ph:
        sum_phi_all = 0
        count = 0
        for ii in range(len(ph[0])):
            if str(ph[col][ii]) != 'nan':
                count+=1
                sum_phi_all = sum_phi_all+np.exp(ph[col][ii]*1j)
        R_all[col] = np.abs(sum_phi_all/count)
    R_all[len(R_all)-1] = R_all[len(R_all)-2]

    phase_coh[kk] = R_all
    
    ph_circ = [((np.arctan2(Y[m, n], X[m, n])+np.pi)/(2*np.pi)) for m in range(Nosc) for n in range(len(t))]
    ph_cell = [((np.arctan2(YY[m, n], XX[m, n])+np.pi)/(2*np.pi)) for m in range(Nosc) for n in range(len(t))]
    hist, xedges, yedges = np.histogram2d(ph_circ, ph_cell, bins=50)
    
    fig = plt.figure(figsize=(12,10))
    plt.imshow(hist.T, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], origin='lower', cmap='plasma', interpolation='bilinear')
    plt.colorbar()
    plt.xlabel(r'Circadian phase ($\theta$/2$\pi$)')
    plt.ylabel(r'Cell cycle phase ($\theta$/2$\pi$)')
    # plt.savefig(output+'Ph_lock_patterns_extr'+str(kappa).replace('.', '')+'_intra'+str(eps).replace('.', '')+'_Nosc'+str(Nosc)+'.svg')
    plt.show()
    
#%%Properties figures -----------------
plt.figure(figsize=(10,8))
for kk in range(len(kappa_vect)):
    plt.plot(t, mean_osc[kk], linewidth=5, label=str(kappa_vect[kk]))    
plt.xlabel('Time [hours]')
plt.xticks([0,24,48,72,96,120])
plt.ylabel('Mean circadian clock oscillation')
plt.legend(loc='best')
# plt.savefig(output+'Mean_phase_circ_changing_extra.svg')
plt.show()

plt.figure(figsize=(10,8))
for kk in range(len(kappa_vect)):
    plt.plot(t, wAn.sinc_smooth(periods[kk], T_c=5), linewidth=5, label=str(kappa_vect[kk]))    
plt.xlabel('Time [hours]')
plt.ylabel('Period')
plt.legend(loc='best')
# plt.savefig(output+'Period_mean_circ_changing_extra.svg')
plt.show()

plt.figure(figsize=(12,10))
for kk in range(len(kappa_vect)):
    plt.plot(t, phase_coh[kk], linewidth=5, label=str(kappa_vect[kk]))
plt.xlabel('Time [hours]')
plt.ylabel('Phase coherence circadian clocks')
plt.legend(loc='best')
# plt.savefig(output+'Ph_coh_diff_extra.svg')
plt.show()



