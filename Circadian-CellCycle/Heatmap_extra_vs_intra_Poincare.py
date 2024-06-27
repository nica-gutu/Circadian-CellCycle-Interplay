#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:06:00 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

def create_custom_colormap(): #input: maximum number divisions
    all_colors = [
    # (220, 20, 60),        # Crimson
    (199, 21, 133),       # Medium Violet Red
    (218, 112, 214),      # Orchid
    (240, 128, 128),      # Light Coral
    (221, 160, 221),      # Plum
    (219, 112, 147),      # Pale Violet Red
    (186, 85, 211),       # Medium Orchid
    (153, 50, 204),       # Dark Orchid
    (102, 51, 153),       # Rebecca Purple
    (147, 112, 219),      # Medium Purple
    (123, 104, 238),      # Medium Slate Blue
    (72, 61, 139)         # Dark Slate Blue
    ]
    
    normalized_colors = [(r / 255, g / 255, b / 255) for r, g, b in all_colors]
    cmap_name = 'custom_colormap'
    
    return LinearSegmentedColormap.from_list(cmap_name, normalized_colors, N=256)

plt.rcParams.update({'font.size': 20})
plt.rcParams['svg.fonttype'] = 'none'

output = '.../'

extra = np.arange(0, 0.041, step=0.001)
intra = np.arange(0, 0.405, step=0.005)
arn_tong = np.zeros((len(extra), len(intra)))
det_xticks = []
coupl_yticks = []

#Coupled oscillators        
gg = 1
ll = 1
a0 = 1
acc = 1
N = 100
dt = 0.2
tf = 200

for mm in range(len(extra)):
    coupl_yticks.append("%.3f"%extra[mm])
    for nn in range(len(intra)):
        mu1, sigma1=20, np.sqrt(4)
        circ_per = np.random.normal(mu1, sigma1, N)
        mu2, sigma2 = 28, np.sqrt(4)
        cell_per = np.random.normal(mu2, sigma2, N)
        det_xticks.append("%.3f"%intra[nn])
        
        kappa = extra[mm]
        eps = intra[nn]
        print(kappa, eps)
        
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
        for m in range(N):
            difference = []
            for n in range(len(t)):
                diff = np.arctan2(Y[m,n], X[m,n])-np.arctan2(YY[m,n], XX[m,n])
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
        R_all[len(t)-1] = R_all[len(t)-2]
        print('mean ph coh', np.mean(R_all[int(len(t)*0.5):len(t)]))
        arn_tong[mm,nn] = np.mean(R_all[int(len(t)*0.5):len(t)])
        
        print('---')
            
arn_tong = np.flipud(arn_tong)

np.savetxt('Arn_tongue_extra_intr2.txt', arn_tong)

#%%Plotting
# file_path = '...'
# arn_tong = np.loadtxt(file_path)
# # arn_tong_dupl = np.hstack((np.fliplr(arn_tong), arn_tong))

# own_colors = ["indigo", "rebeccapurple", "mediumorchid", "orchid", "violet", "lightpink", "pink", "lavenderblush"]
# cmap = mcolors.LinearSegmentedColormap.from_list("", own_colors)

# fig = plt.figure(figsize=(12, 8))
# plt.xticks([0, 20, 40],[0, 0.1, 0.2])
# color_map = plt.imshow(arn_tong,  interpolation='bicubic', aspect='auto')#, norm=PowerNorm(gamma=5))
# color_map.set_cmap(cmap)
# plt.xlabel('Intracellular coupling')
# plt.yticks([0, 10, 20, 30, 40, 50],[0.05, 0.04, 0.03, 0.02, 0.01, 0])
# plt.ylabel('Extracellular coupling')
# ax = plt.gca()
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.1)
# cb = plt.colorbar(cax=cax)
# cb.set_label('Mean phase coherence of phase difference')
# plt.savefig(output+'Arn_tongue_extra_intra.svg')
# plt.show()


