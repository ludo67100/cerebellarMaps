# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:01:30 2020

@author: Ludovic.SPAETH
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sn

import pandas as pd 
import numpy as np

from matplotlib import pyplot as plt


url = 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max'
savedir= 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max/FIGURES'

conditions = ['WT','ENR','LC','LS','EC','ES']
colors = ['skyblue','green','lightcoral','black','orange','purple']
cmaps = ['Blues','Greens','Reds','Greys','Oranges','Purples']

winsorize = False
cutoff = 0.05 #For winsorize

N = 8
r_min = 0.55 #Minimum r value for detection as correlation,
zscore = 3.09

savefig = False
savedata = False

vmin = r_min

individual_plot = False 

fig3, ax = plt.subplots(2,1,figsize=(N,N))
ax[1].set_xticks(range(len(conditions)))
ax[1].set_xticklabels(conditions)
ax[1].set_ylabel('Count')

#-------------------------------------THE SERIOUS SHIT-----------------------------------

COUNTER = []

for cond,color,cmap,idx in zip(conditions,colors,cmaps,range(len(conditions))) :
    
    print ('--- Condition:{} ---'.format(cond))

    Rs = np.zeros((N,N))
    
    fig2,axx = plt.subplots(1,1,figsize=(N,N))
    
    path = '{}/{}_Map_2D_Average_Amp_per_zones_wNANS_{}.xlsx'.format(url,cond,zscore) 

    #Get dataset, removing first column with names
    dataset = pd.read_excel(path).iloc[:,1:]

    bands = dataset.columns
    animals = dataset.index
    
    #Do pearson correlation
    corrMatrix = dataset.corr(method='pearson')
    
    #Count how much time a given zone is correlated to another
    countMatrix = corrMatrix[corrMatrix.abs()>=r_min].count(axis=0)-1
    COUNTER.append(countMatrix)
    
    axx=sn.heatmap(corrMatrix,annot=True,cmap=cmap,vmin=r_min,vmax=0.9,
                   cbar_kws={'label': 'Pearson coeff'}).set_title(cond)
    
    list_of_pearsons_coeff = [x for x in np.unique(corrMatrix) if x !=1.]
    
    ax[0].hist(list_of_pearsons_coeff,bins=10,cumulative=True,density=True,label=cond, color=color,histtype='step',linewidth=2)
    ax[0].set_xlabel('Pearson coeff')
    ax[0].set_ylabel('Normalized count')
    ax[0].legend(loc='best')

    
    #Determine how many zones are correlated within each condition 
    positive_count = 0
    negative_count = 0    

    for coeff in list_of_pearsons_coeff :
        if coeff >= r_min:
            positive_count += 1
            
        elif coeff <= -r_min:
            negative_count -=1
            
        else:
            continue
            
    print('Postive and negative counts = {}/{}'.format(positive_count,negative_count))
            
    
    ax[1].bar(idx,positive_count,color=color,label=cond)
    ax[1].bar(idx,negative_count,color=color,alpha=0.5)
    
    #3D surface plot---------------------------------------------------------- 
    figure3D = plt.figure(figsize=(8,6))
    ax3D = plt.subplot(111,projection='3d')
    
#    # Remove gray panes and axis grid
#    ax3D.xaxis.pane.fill = False
#    ax3D.xaxis.pane.set_edgecolor('white')
#    ax3D.yaxis.pane.fill = False
#    ax3D.yaxis.pane.set_edgecolor('white')
#    ax3D.zaxis.pane.fill = False
#    ax3D.zaxis.pane.set_edgecolor('white')
#    ax3D.grid(False)
#    
#    # Remove z-axis
#    ax3D.w_zaxis.line.set_lw(0.)
#    ax3D.set_zticks([])
    
    # Create meshgrid
    X, Y = np.meshgrid(np.linspace(0, 2, len(corrMatrix)), np.linspace(0, 2, len(corrMatrix)))
    
    #Smooooooth it !
    xnew, ynew = np.mgrid[-1:1:20j, -1:1:20j]
    
    from scipy import interpolate
    tck = interpolate.bisplrep(X, Y, corrMatrix, s=0)
    znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    
    # Plot surface
    plot = ax3D.plot_surface(X=xnew, Y=ynew, Z=znew, cmap='bwr', vmin=-1, vmax=1,antialiased=True)
    ax3D.view_init(azim=45,elev=45)
    
    if savefig==True:
        fig2.savefig('{}/{}_Pearson_Matrix_zscore{}.pdf'.format(savedir,cond,zscore))
        fig2.savefig('{}/{}_Pearson_Matrix_zscore{}.png'.format(savedir,cond,zscore))  
        
ax[1].plot(range(len(conditions)),np.zeros(len(conditions)),color='black',linestyle='--')

if savefig==True:
    fig3.savefig('{}/Group_Data_zscore{}.pdf'.format(savedir,zscore))
    fig3.savefig('{}/Group_Data_zscore{}.png'.format(savedir,zscore))
    
    
    
#Counting correlations 
    
countDf = pd.concat(COUNTER,axis=1)
countDf.columns=conditions


finalFig, axis = plt.subplots(1,2)
plt.suptitle(url.split('/')[-1])
axis[0].set_title('Sub-zones')
axis[0].set_ylabel('# of corr. with other zones (R>{})'.format(r_min))
axis[1].set_title('Zones')


conditions = ['WT','ENR','LC','LS','EC','ES']
colors = ['skyblue','green','lightcoral','black','orange','purple']

countDf.plot.bar(ax=axis[0], color=colors)


#Merge zones to get A, AX and B zones only

Bzone =  countDf.T['B_contra'] + countDf.T['B_ipsi']
AXzone = countDf.T['Ax_Contra'] + countDf.T['Ax_ipsi']
Azone = countDf.T['A_lat_contra'] + countDf.T['A_med_contra'] + countDf.T['A_med_ipsi'] + countDf.T['A_lat_ipsi']

moduleDf = pd.concat([Bzone, AXzone, Azone], axis=1)
moduleDf.columns=['B','AX','A']

moduleDf.T.plot.bar(ax=axis[1],colors=colors)

plt.tight_layout()





    
    
    
    
    
    
        
        
