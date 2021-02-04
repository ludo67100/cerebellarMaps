# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 16:48:37 2020

@author: ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd 
import numpy as np 


#The file containing linear reg p values for each zone 
file = 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max/FIGURES/00_SORTED_PVALUES_MICROZONES.xlsx' 

zones = ['B_contra','Ax_Contra','A_lat_contra','A_med_contra',
         'A_med_ipsi','A_lat_ipsi','Ax_ipsi','B_ipsi']

saveDir = 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max/FIGURES'

colors = ['skyblue','limegreen','green','lightcoral','black','orange','purple']

df = pd.read_excel(file,header=0,index_col=0)


globalDf = []

for i in range(df.shape[0]): 
    
    dataCond = df.iloc[i,:]
    
    temp = []
    
    for index in range(len(dataCond.index)): 
        
        data = df.iloc[i,index]
        
        if data <= 0.05: 
            temp.append(1)
            
        else:
            temp.append(0)
            
    pvalues = pd.DataFrame(temp).values.reshape(8,8)
    
    sumOfPvalues = pd.DataFrame(np.nansum(pvalues,axis=1)-1,index=zones)
    
    globalDf.append(sumOfPvalues)
    


counts = pd.concat(globalDf, axis=1)
counts.columns=df.index


counts.to_excel('{}/01_Microzones_Connectivity.xlsx'.format(saveDir))

#Do a heatmap 
from matplotlib import pyplot as plt
import seaborn as sn 

fig, ax = plt.subplots(1,1,figsize=(3,6))
sn.heatmap(counts.T)

fig.savefig('{}/Zone_count_heatmap.pdf'.format(saveDir))

#Add somehow of histgoram-ridge plot
fig, ax = plt.subplots(counts.shape[1],1,sharex=True,sharey=True,figsize=(3,6))

for i in range(counts.shape[1]): 
    
    ax[i].grid(which='both',axis='y')
    ax[i].bar(np.arange(0,len(counts.iloc[:,i]),1),counts.iloc[:,i],color=colors[i],label=counts.iloc[:,i].name)
    ax[i].set_ylabel('count')
    ax[i].legend(loc='best')

    
    if i == counts.shape[1]-1:
        
        ax[i].set_xticks(np.arange(0,len(zones),1))
        ax[i].set_xticklabels(zones,rotation='45')
        
fig.savefig('{}/Zone_count_barplots.pdf'.format(saveDir))
    
    





            