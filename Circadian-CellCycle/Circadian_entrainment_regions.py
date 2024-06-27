#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:30:00 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 26})
plt.rcParams['svg.fonttype'] = 'none'

output = '.../'
path2 = '.../'

dt = 0.1 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

coupling = np.arange(0, 0.0505, step=0.0005)
dispersion = np.arange(0, 10.1, step=0.1)
arn_tong = np.zeros((len(coupling),len(dispersion)))

#Simulation parameters        
N = 20
dt = 0.1
tf = 150

#Coupled oscillators parameters
gg = 1
ll = 1
a0 = 1
acc = 1

eps = 0 #Intracellular coupling

for mm in range(len(coupling)):
    for nn in range(len(dispersion)):
        delta = dispersion[nn]
        print('dispersion', delta)

        mu1, sigma1 = 24, np.sqrt(delta)
        circ_per=np.random.normal(mu1, sigma1, N)
                
        kappa = coupling[mm]
        print(kappa)
        
        t = np.arange(0, tf, step=dt)
        X = np.zeros((N, len(t)))
        Y = np.zeros((N, len(t)))
        
        X[:,0] = np.random.uniform(-1, 1, size=(N))
        Y[:,0] = np.random.uniform(-1, 1, size=(N))
            
        for i in range(len(t)-1):
            sum_x = 0
            sum_y = 0
            for m in range(N):                
                middle = -gg*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0)
                X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/circ_per[m]+sum_x*kappa/(2*N))
                Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/circ_per[m]+sum_y*kappa/(2*N))
                
                for k in range(N):
                    sum_x = sum_x+X[k,i] #global coupling
                    sum_y = sum_y+Y[k,i] #global coupling

#%%Synchronization index
        mean_cell_sqr_time = np.mean((np.mean(X, axis=0))**2)
        mean_cell_time_sqr = (np.mean(np.mean(X, axis=0)))**2
        avg = 0
        for i in range(N):
            avg = avg+np.mean((X[i,:])**2)-(np.mean(X[i,:]))**2
        R = (mean_cell_sqr_time-mean_cell_time_sqr)/(avg/N)
        arn_tong[mm,nn] = R #np.mean(R_all[int(len(t)*0.75):len(t)]) #to save
        
            
arn_tong = np.flipud(arn_tong)

np.savetxt(path2+'Arn_tongue_extracell_coupling.txt', arn_tong)

file_path = path2+'Arn_tongue_extracell_N20_k005.txt'
arn_tong = np.loadtxt(file_path)
arn_tong_dupl = np.hstack((np.fliplr(arn_tong), arn_tong))

fig = plt.figure(figsize=(14, 10))
plt.xticks(np.linspace(0, 200, 9), [10, 7.5, 5, 2.5, 0, 2.5, 5, 7.5, 10])
color_map = plt.imshow(arn_tong_dupl, interpolation='bicubic', aspect='auto')
color_map.set_cmap('jet')
plt.xlabel('Dispersion circadian periods [hours]')
plt.yticks(np.linspace(0,100,11), np.round(np.linspace(0.05, 0, 11), 5))
plt.locator_params(axis='y', nbins=20)
plt.ylabel('Extracellular coupling')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cb = plt.colorbar(cax=cax)
cb.set_label('Synchronization index R')
# plt.savefig(output+'Extra_coupled_circadian_arntongue.svg')
plt.show()
