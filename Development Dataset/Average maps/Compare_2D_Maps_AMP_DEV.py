# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 12:16:15 2020

@author: Ludovic.SPAETH
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
import scipy.ndimage as ndim
import seaborn as sn
import pandas as pd 
import math

#Where are the data ?
datasource = 'E:/000_PAPER/Development/Amplitude_Analysis/02_AVERAGE_2D'

#Where do we save the data ?
savedir = 'E:/000_PAPER/Development/Amplitude_Analysis/02_AVERAGE_2D'

#What kind of data do we use ?
#data = '2D_MedianZscore.csv'
data = '2D_AvgAmplitude.csv'

#Zebrin file 
zebrinFile = 'E:/000_PAPER/Mesures_ZII_LowRes_Adult_and_Dev.xlsx'

#Define the groups
groupsToCompute = 4
groups = ['P9P10','P12P13','P14P18','P30P40']
colors = ['lightskyblue','skyblue','deepskyblue','royalblue']

#Pairs to analyse for stats
pairs= [("P9P10", "P12P13"), ("P9P10", "P14P18"), ("P9P10", "P30P40"), ("P12P13", "P14P18"),
        ("P12P13", "P30P40"),("P14P18", "P30P40")]

#Half width of gaussian convolution kernel 
sigma = 18

SDfactor = 3  #3 normalement

#Vmin vmax for 2D maps 
#Vmin is threshold based on avg noise value in Noise_Values_Avg_and_Median.xlsx file
# vmin = avg_noise + X*noise SD
vmin, vmax = 4.38+SDfactor*5.01,60

#Vmin and Vmax for differntial plot
delta_vmin, delta_vmax = -30,30

#Save the figures ? 
saveFig = False

#Save the data ? 
saveData = False

#Store the convolved maps
convolvedMaps = []

#Functions---------------------------------------------------------------------
def count_event_mediolateral_axis(map2D, threshold):
    '''
    Returns hist of events (if site > threshold) along the mediolateral axis 
    '''
    events = []
    
    for i in range(map2D.shape[1]): 
    
        col = map2D[:,i]
        
        if math.isnan(col[0]) == True:    
            count = 0
        
        else:
            count = [1. for x in col if x > threshold]
        
        events.append(np.nansum(count)/map2D.shape[0]*100)
        
    return np.asarray(events)


def count_event_anteroposterior_axis(map2D, threshold):
    '''
    Returns hist of events (if site > threshold) along the anteroposterios axis 
    '''
    events = []
    
    for i in range(map2D.shape[0]): 
    
        col = map2D[i,:]
        
        if math.isnan(col[0]) == True:    
            count = 0
        
        else:
            count = [1. for x in col if x > threshold]
        
        events.append(np.nansum(count)/map2D.shape[1]*100)
        
    return np.asarray(events)
    

#Create kernel for convolution
# First a 1-D  Gaussian
t = np.linspace(-10, 10, sigma)
bump = np.exp(-0.5*t**2)
bump /= np.trapz(bump) # normalize the integral to 1

# make a 2-D kernel out of it
kernel = bump[:, np.newaxis] * bump[np.newaxis, :]

flattened_maps = pd.DataFrame()

#Iterate on each group to compute convolution
for group, color in zip(groups,colors):
    print ('{} 2D map : convolved'.format(group))
    
    datadir = '{}/{}'.format(datasource,group)
    
    zebrinDf = pd.read_excel(zebrinFile,header=0,index_col=0,sheet_name=group)
    
    zebrinValues = zebrinDf.loc['MEAN (normalized)'].values
    
    list_of_2D_files = ['{}/{}'.format(datadir,x) for x in os.listdir(datadir) if x.endswith('{}'.format(data))]
     
    #Get the raw map
    rawMap = np.abs(np.genfromtxt('{}/{}_{}'.format(datadir,group,data),
                            delimiter=','))
    
    
    np.savetxt('{}/{}_rawmap.csv'.format(savedir,group),rawMap,delimiter=',')
    
    
    #Get positional array
    positions = [int(x) for x in np.genfromtxt('{}/{}_POSITIONAL_ARRAY.csv'.format(datadir,group))]
    
#    #Convolute map with gaussian kernel
#    convolvedMap = ndim.gaussian_filter(rawMap,
#                                        sigma=sigma,
#                                        order=0,
#                                        mode='constant')
    
    #Convolute map with gaussian kernel
    convolvedMap = ndim.convolve(rawMap,
                                 weights=kernel,
                                 mode='constant',cval=vmin)
    
    flattened_maps[group] = np.ravel(rawMap)
    
    #Compute histogram on both axis of the map
    mediolateralAxis = count_event_mediolateral_axis(convolvedMap,threshold=vmin)
    anteropostAxis = count_event_anteroposterior_axis(convolvedMap,threshold=vmin)
    
    
    #Turn map in DataFrame
    df = pd.DataFrame(convolvedMap, index=['40','80','120','160'], columns=positions[:-1])
    
    if saveData == True:
        
        df.to_excel('{}/{}_convolved_2D_map.xlsx'.format(savedir, group))
    
    #Append for later
    convolvedMaps.append(convolvedMap)    
    
    #Plot plot plot
    #First the maps
    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(10,4))
    sn.heatmap(rawMap, ax=ax1,cbar_kws={'label': 'Syn. Amp.'}, vmin=vmin, vmax=vmax,
               xticklabels=positions).set_title('{} Raw Map'.format(group))
    plt.xticks(rotation=90)  
    
    sn.heatmap(convolvedMap, ax=ax2,cbar_kws={'label': 'conv. Amp'}, vmin=vmin, vmax=vmax,
               xticklabels=positions).set_title('{} Convolved Map (Gaussian kernel, sigma={})'.format(group,sigma))
    
    
    plt.xticks(rotation=90)   
    plt.tight_layout(h_pad=2,w_pad=None)
    
    #The the histograms 
    figg, hist = plt.subplots(1,2,figsize=(16,4))
    figg.suptitle('{} axis distribution (threshold={}*noise SD)'.format(group,SDfactor))

    hist[0].set_ylabel('% in column')
    hist[0].set_ylim(0,100)
    hist[0].set_xlabel('Mediolateral axis #')
    hist[0].bar(np.arange(0,len(mediolateralAxis),1),mediolateralAxis,color=color)  
    
    hist[0].set_xticks(np.arange(0,len(mediolateralAxis)))
    hist[0].set_xticklabels(positions, rotation='90')
    
    #Add zebrins on mediolateral histogram plot
    #hist[0].axvspan(zebrinValues[1],zebrinValues[2],color='green',alpha=0.2)
    
    
    hist[1].set_ylabel('AnteroPosterior axis #')
    hist[1].set_xlabel('% in row')
    hist[1].set_xlim(0,100)
    

    
    depths = np.flip(['0-40','40-80','80-120','120-160'],axis=0)
    
    hist[1].barh(np.arange(0,len(anteropostAxis),1),np.flip(anteropostAxis,axis=0),color=color)   

    hist[1].set_yticks(np.arange(0,len(anteropostAxis),1))
    hist[1].set_yticklabels(depths)
    
    
    if saveFig == True:
    
        fig.savefig('{}/{}_convolved_2D_map.pdf'.format(savedir,group))
        fig.savefig('{}/{}_convolved_2D_map.png'.format(savedir,group))

        figg.savefig('{}/{}_axis_histograms.pdf'.format(savedir,group))
        figg.savefig('{}/{}_axis_histograms.png'.format(savedir,group))

#Substract WT pattern to other conditions
fig2,axx = plt.subplots(len(groups),1,sharex=True,sharey=True,figsize=(10,10))
for i in range(len(groups)):
    
    datadir = '{}/{}'.format(datasource,groups[i])
    
    #Get positional array
    positions = [int(x) for x in np.genfromtxt('{}/{}_POSITIONAL_ARRAY.csv'.format(datadir,groups[i]))]

    
    wtConvolvedMap = convolvedMaps[0]
    
    differentialMap = ndim.convolve(convolvedMaps[i]-wtConvolvedMap, weights=kernel, mode='constant')
    
    chart = sn.heatmap(differentialMap,ax=axx[i],annot=False, vmin=delta_vmin, vmax=delta_vmax,
                       cbar_kws={'label': 'Delta'},cmap='coolwarm',
                       xticklabels=positions).set_title('{} - WT'.format(groups[i]))

plt.xticks(rotation=90)        
plt.tight_layout(h_pad=2,w_pad=None)


#TODO : save stuff
if saveFig == True:
    
    fig2.savefig('{}/Comparison.pdf'.format(savedir))
    fig2.savefig('{}/Comparison.png'.format(savedir))    
    
    
import pingouin as pg 
from scipy import stats
    
#Now stats and group data
mainfig, mainplot = plt.subplots(1,1)
mainfig.suptitle('Average Amplitudes distributions')
mainplot.set_ylabel('Average Amplitudes (pA)')

sn.boxplot(data=flattened_maps,ax=mainplot,palette=colors)
sn.swarmplot(data=flattened_maps,ax=mainplot,color='black',size=2) 

normality = pg.normality(flattened_maps)

if normality['normal'][0] == False:
    groupStat = stats.kruskal(flattened_maps.values[:,0],
                              flattened_maps.values[:,1],
                              flattened_maps.values[:,2],
                              flattened_maps.values[:,3])
    print ('Normality failed, KW test (pvalue)={}'.format(groupStat[1]))
    
    if groupStat[1] < 0.05:
        
        for group in groups:
            
            controlGroup = flattened_maps['P30P40'].values
            compareGroup = flattened_maps[group].values
            
            mwu_test = stats.mannwhitneyu(controlGroup,compareGroup,alternative='two-sided')
            print ('MWU test: P30P40 vs {} (corrected p-value)= {}'.format(group,mwu_test[1]))
            
            if mwu_test[1] < 0.05/groupsToCompute-1:
                print('Difference is statistical (correction applied)')
            else:
                print ('No statistical difference after correction')

else:
    groupStat = stats.anova(flattened_maps.values[:,0],
                              flattened_maps.values[:,1],
                              flattened_maps.values[:,2],
                              flattened_maps.values[:,3])
    print ('Normality is confirmed, Anova (pvalue)={}'.format(groupStat[1]))

    


