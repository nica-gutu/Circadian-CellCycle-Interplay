#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:10:05 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

output = '.../'

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

kappa_vect = [0, 0.05] #extracellular circadian coupling
eps_vect = [0, 0.15] #intracellular circadian - cell cycle #0.1 or 0.15

plt.figure(figsize=(12,10))

for kappa in kappa_vect:
    for eps in eps_vect:
        #Vectors initialization
        t = np.arange(0,tf,step=dt)
        X = np.zeros((Nosc,len(t)))
        XX = np.zeros((Nosc,len(t)))
        Y = np.zeros((Nosc,len(t)))
        YY = np.zeros((Nosc,len(t)))
        
        #Initial position
        # np.random.seed(0)
        X[:,0] = np.random.uniform(-1,1,size=(Nosc))
        XX[:,0] = np.random.uniform(-1,1,size=(Nosc)) 
        Y[:,0] = np.random.uniform(-1,1,size=(Nosc))
        YY[:,0] = np.random.uniform(-1,1,size=(Nosc))
        
        #Period distribution
        mu1, sigma1 = 22, np.sqrt(6)
        period_circ = np.random.normal(mu1, sigma1, Nosc)
        mu2, sigma2 = 26, np.sqrt(6)
        period_cell = np.random.normal(mu2, sigma2, Nosc)
        
        
        #Integration over time and population
        for i in range(len(t) - 1):
            sum_x = 0
            sum_y = 0
            
            for m in range(Nosc):
                middle = -gg*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0)
                X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/period_circ[m]+sum_x*kappa/(2*Nosc))
                Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/period_circ[m]+sum_y*kappa/(2*Nosc))
                
                mid = -ll*(np.sqrt(XX[m,i]**2+YY[m,i]**2)-acc)
                XX[m,i+1] = XX[m,i]+dt*(mid*XX[m,i]-2*np.pi*YY[m,i]/period_cell[m]+eps*(X[m,i]+XX[m,i])/2)
                YY[m,i+1] = YY[m,i]+dt*(mid*YY[m,i]+2*np.pi*XX[m,i]/period_cell[m]+eps*(Y[m,i]+YY[m,i])/2)
        
                for k in range(Nosc):
                    sum_x = sum_x+X[k,i] 
                    sum_y = sum_y+Y[k,i] 
                    
        #Phase coherence of phase differences    
        ph_diff = pd.DataFrame(columns=np.arange(0, Nosc, step=1))
        for m in range(Nosc):
            difference=[]
            for n in range(len(t)):
                diff = np.arctan2(Y[m,n], X[m,n])-np.arctan2(YY[m,n], XX[m,n]) #which oscillator
                phasediff = np.arctan2(np.sin(diff), np.cos(diff))
                difference.append(phasediff)   
            ph_diff[m] = difference
              
        ph_diff = ph_diff.T 
        R_all = np.zeros(len(t))
        for col in ph_diff:
            sum_phi_all = 0
            count = 0
            for ii in range(len(ph_diff[0])):
                if str(ph_diff[col][ii]) != 'nan':
                    count+=1
                    sum_phi_all = sum_phi_all+np.exp(ph_diff[col][ii]*1j)
            R_all[col] = np.abs(sum_phi_all/count)
        R_all[len(R_all)-1] = R_all[len(R_all)-2]
        
        plt.plot(t, R_all, linewidth=5, label='k:'+str(kappa)+' eps:'+str(eps))
plt.xlabel('Time [hours])')
plt.ylabel('Phase coherence of phase differences')
plt.legend(loc='best')
plt.savefig(output+'Ph_coh_ph_diff_4extremecases.svg')
plt.show()

