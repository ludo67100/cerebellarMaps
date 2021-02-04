# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:57:42 2020

@author: Ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42



import pandas as pd 
import numpy as np 
    


vmin = 4.38+2*5.01


conditions = ['P9P10','P12P13','P14P18','P30P40']

from matplotlib import pyplot as plt 

fig, ax = plt.subplots(4,1)

for cond,i in zip(conditions,range(len(conditions))): 
    
    map2d = pd.read_excel('E:/000_PAPER/Development/Amplitude_Analysis/02_AVERAGE_2D/{}_convolved_2D_map.xlsx'.format(cond))
    rawmap = np.genfromtxt('E:/000_PAPER/Development/Amplitude_Analysis/02_AVERAGE_2D/{}_rawmap.csv'.format(cond),delimiter=',')


    plot = ax[i].imshow(rawmap, interpolation='sinc', cmap='magma',vmin=vmin,vmax=75)
    ax[i].set_title(cond)
    
    fig.colorbar(plot,ax=ax[i])
    
    
    
# for figure 2  in paper, put individual maps in magma cmap and scale color bars
    

files = ['E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS/P9P10/150410(2)/150410(2)_Amp_2D_OK.csv',
         'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS/P12P13/150429(1)/150429(1)_Amp_2D_OK.csv',
         'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS/P14P18/151117(1)/151117(1)_Amp_2D_OK.csv',
         'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS/P30P40/150423(1)/150423(1)_Amp_2D_OK.csv']

fig, ax = plt.subplots(4,1)

for file, i in zip(files, range(len(files))):

    rawmap = np.abs(np.genfromtxt(file,delimiter=','))
    
    plot = ax[i].imshow(rawmap, interpolation=None, cmap='magma',vmax=100)
    ax[i].set_title(file.split('/')[-1])
    
    fig.colorbar(plot,ax=ax[i])

    
    
    
