#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:55:30 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def create_custom_colormap(num_colors): #input: maximum number divisions
    all_colors = [
        (135/255, 206/255, 250/255),  # skyblue
        # (255/255, 228/255, 181/255),  # moccasin
        (255/255, 228/255, 225/255),     #mistyrose
        (216/255, 191/255, 216/255),  # thistle
        (240/255, 128/255, 128/255),  # lightcoral
        (95/255, 158/255, 160/255),   # cadetblue
        (255/255, 215/255, 0/255),    # gold
    ]
    
    if num_colors not in range(3, 7):
        raise ValueError("Number of colors must be between 3 and 6 inclusive")

    colors = all_colors[:num_colors]
    
    cmap_name = 'custom_colormap'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=256)


path = '/Users/nicagutu/Nextcloud/Manuscripts/Circadian-CellCycle/Data/1stExp_ilastik/Division_matrix/'
output = '/Users/nicagutu/Nextcloud/Manuscripts/Circadian-CellCycle/Figures/Fig4/'

dose = ['untreated', '5uM', '10uM'] 
density = ['high', 'medium', 'low'] 
channel2 = 'cell_cycle'

dt = 0.5 # the sampling interval, 0.5hours
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

for i in density:
    for j in dose:

        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Division matrix
        data = pd.read_excel(path+'divisions_'+str(condition)+r'.xlsx')
        
        IDs = data.columns
        time_points = len(data.index.values) 
        properties = {'num_divisions':[],'time_1st':[]}
        
        for col in IDs:
            properties['num_divisions'].append(sum(data[col]))
            count=0
            for jj in range(len(data[col])):
                if sum(data[col]) == 0 and count <1:
                    count+=1
                    properties['time_1st'].append(0)
                elif sum(data[col]) != 0:
                    if data[col][jj] == 1 and count <1:
                        properties['time_1st'].append(data.index.values[jj])
                        count+=1
                    
        properties['num_divisions'],properties['time_1st'],IDs=zip(*sorted(zip(properties['num_divisions'],properties['time_1st'],IDs)))
        
        divisions = data.reindex(columns=IDs)
        division_profile = (divisions.to_numpy())
        division_profile = division_profile.T

        for xx in range(len(IDs)):
            divisions = []
            for yy in range(time_points):
                if division_profile[xx, yy] == 1:
                    divisions.append(yy)
        
            if len(divisions) > 0:
                for n in range(len(divisions) - 1):
                    division_profile[xx, divisions[n]+1:divisions[n+1]] = n+1
        
                division_profile[xx, divisions[-1]+1:time_points] = len(divisions)
        
                if len(divisions) > 1:
                    division_profile[xx, divisions[-1]] = len(divisions)
                    
        plt.figure(figsize=(10, 20))
        positions = (0, 100, 200)
        label_pos = (0, 50, 100)
        plt.xticks(positions, label_pos)
        color_map = plt.imshow(division_profile)
        color_map.set_cmap(create_custom_colormap(max(properties['num_divisions'])+1))
        plt.xlabel('Time(h)')
        plt.ylabel('Cell number')
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cb = plt.colorbar(cax=cax)
        labels_list = np.arange(0, sum(divisions), 1)
        loc = labels_list+0
        cb.set_ticks(loc)
        cb.set_ticklabels(labels_list)
        cb.set_label('No. divisions', rotation=270, labelpad=25)
        plt.savefig(output+'Division_profile_'+str(condition)+'.svg', format='svg')
        plt.show()
        
        