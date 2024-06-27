#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:49:54 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 20})
plt.rcParams['svg.fonttype'] = 'none'

output = '.../'

coupling = np.arange(0, 0.42, step=0.02)
detuning = np.arange(0, 8.1, step=0.1)

arn_tong = np.zeros((len(coupling), len(detuning)))
det_xticks = []
coupl_yticks = []

#Coupled oscillators parameters
N = 200
dt = 0.1
tf = 200

gg = 1
ll = 1
a0 = 1
acc = 1

kappa = 0.1 #extracellular coupling

for mm in range(len(coupling)):
    for nn in range(len(detuning)):
        delta = detuning[nn]
        mu1, sigma1 = 20+delta, np.sqrt(4)
        circ_per = np.random.normal(mu1, sigma1, N)
        mu2, sigma2 = 28-delta, np.sqrt(4)
        cell_per = np.random.normal(mu2, sigma2, N)

        overlap = NormalDist(mu=mu1, sigma=sigma1).overlap(NormalDist(mu=mu2, sigma=sigma2))
        print('overlap', overlap*100)
        
        eps = coupling[mm]
        print('coupling, detuning:', eps, mu2-mu1)

        t = np.arange(0, tf, step=dt)
        X = np.zeros((N, len(t)))
        Y = np.zeros((N, len(t)))
        XX = np.zeros((N, len(t)))
        YY = np.zeros((N, len(t)))
        
        X[:,0] = np.random.uniform(-1, 1, size=(N))
        Y[:,0] = np.random.uniform(-1, 1, size=(N))
        XX[:,0] = np.random.uniform(-1, 1, size=(N))
        YY[:,0] = np.random.uniform(-1, 1, size=(N))
            
        for i in range(len(t)-1):
            sum_x = 0
            sum_y = 0
            for m in range(N):                
                middle = -gg*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0)
                X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/circ_per[m]+sum_x*kappa/(2*N))
                Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/circ_per[m]+sum_y*kappa/(2*N))
                
                mid = -ll*(np.sqrt(XX[m,i]**2+YY[m,i]**2)-acc)
                XX[m,i+1] = XX[m,i]+dt*(mid*XX[m,i]-2*np.pi*YY[m,i]/cell_per[m]+eps*(X[m,i]+XX[m,i])/2)
                YY[m,i+1] = YY[m,i]+dt*(mid*YY[m,i]+2*np.pi*XX[m,i]/cell_per[m]+eps*(Y[m,i]+YY[m,i])/2)        
                for k in range(N):
                    sum_x = sum_x+X[k,i] #global coupling
                    sum_y = sum_y+Y[k,i] #global coupling
                
        ph_diff = pd.DataFrame(columns=np.arange(0, N, step=1))
        std_all = []
        count_sync = 0
        for m in range(N):
            difference = []
            for n in range(len(t)):
                diff = np.arctan2(Y[m,n], X[m,n])-np.arctan2(YY[m,n], XX[m,n])
                phasediff = np.arctan2(np.sin(diff), np.cos(diff))
                difference.append(phasediff)   
            ph_diff[m] = difference
            std_ph = np.std(difference[int(len(t)/2):len(t)])
            std_all.append(std_ph)
            if std_ph < 0.5:
                count_sync+=1
        arn_tong[mm,nn] = count_sync
        
        print('synch',count_sync)
        print('---')
        
arn_tong = np.flipud(arn_tong)
np.savetxt('Arn_tongue_intr_N200_k005.txt', arn_tong)

file_path = '/Users/nicagutu/Nextcloud/Manuscripts/Circadian-CellCycle/Code/Poincare_model/Arn_tongue_intr_N200_k005.txt'
arn_tong = np.loadtxt(file_path)

fig = plt.figure(figsize=(12, 10))
color_map = plt.imshow(arn_tong, interpolation='bicubic', aspect='auto')
color_map.set_cmap('jet')
plt.yticks([0, 10, 20], [0.4, 0.2, 0])
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80], [-8, -6, -4, -2, 0, 2, 4, 6, 8])
plt.xlabel('Detuning')
plt.ylabel('Intracellular coupling')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cb = plt.colorbar(cax=cax)
cb.set_label('# locked oscillators')
# plt.savefig(output+'Arn_tongue_intra'+str(k)+'.svg')
plt.show()
