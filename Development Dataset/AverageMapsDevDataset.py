# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:22:58 2020

This code computes median synaptic patterns and median maps from individual synaptic maps
AMPLITUDE way : inputs are amplitude (pA) based maps 

@author: ludov
"""

import numpy as np 
from numpy import genfromtxt as gen
from matplotlib import pyplot as plt 
import os 
import pandas as pd 
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

#The groups to analyse
conditions = ['P9P10','P12P13','P14P18','P30P40']
colors = ['lightskyblue','skyblue','deepskyblue','royalblue']

#The sheetnames in excel file
sheets=['P9P10','P12P13','P14P18','P30P40']


#General directory to find the data
inputDir = 'D:/03_FORMATED_DATA/For Paper/EPHYS/Development_Dataset'

#Where to save datas/figures
outputDir = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/REVISION CODE/00_FINAL_CODE/Figure2/2D/Average Maps'

#Zebrin file
zebrinFile = 'D:/03_FORMATED_DATA/For Paper/EPHYS/Mesures_ZII_LowRes_Adult_and_Dev.xlsx'

#To constrict 1D plots
ylim = 20

binForMedian = 20  #In % of P1- : 10 is the last value used

left,right = 210,210 #Borders in %P1- 

zscoreLimit =2.0 #Limit of significance for zscore

vmin, vmax = -150,-10 #For maps plot to align every conditions
minz, maxz = 3.09,20

#Interpolation ?
interpolationType = 'sinc'

#Do we save anything ?
saveFig = True #For figures
saveData = True


#------------------------------------FUNCTIONS---------------------------------------------------------
#---------------------------------DO NOT MODIFY--------------------------------------------------------
def MAD(a,axis=None):
    '''
    Computes median absolute deviation of an array along given axis
    '''
    #Median along given axis but keep reduced axis so that result can still broadcast along a 
    
    med = np.nanmedian(a, axis=axis, keepdims=True)
    mad = np.median(np.abs(a-med),axis=axis) #MAD along the given axis
    
    return mad 

def stack_lines(_list):
    
    stacked_map = np.vstack((_list[0],
                            _list[1],
                            _list[2],
                            _list[3]))
    
    return stacked_map


#Iterate for each condition
for condition,sheet,i in zip(conditions,sheets,range(len(conditions))): 
    
    parentDir = '{}/{}'.format(inputDir,condition)
    print(parentDir)
    
    listOfExperiments = os.listdir(parentDir)
    
    #Load zebrin file (i.e. corresponding datasheet)
    zebrins = pd.read_excel(zebrinFile,sheet_name=sheet,
                            index_col=0,header=1)
    
    #Matrix to append maps and positions
    H = 4           #Map height in sites
    L = 32          #Map width in sites
    N = len(listOfExperiments)
    
    _mat = np.zeros((H,L,N,3)) #[0] for map, [1] for position and [2] for Zscore


    #Iterate on each experiment
    for experiment,idx in zip(listOfExperiments,range(N)):
        print (experiment)        
        manipPath = '{}/{}'.format(parentDir,experiment)
        
        #Load map in Matrix
        _mat[:,:,idx,0]=gen('{}/{}_Amp_2D_OK.csv'.format(manipPath,experiment),delimiter=',')
        
        #Get the positions
        pos = gen('{}/{}_Positions_cp_centered_OK.csv'.format(manipPath,experiment),delimiter=',')                
        pos_2D = (pos,pos,pos,pos)        
        _mat[:,:,idx,1]=np.reshape(pos_2D,(H,L))
        
        #And now the 2D Zscore    
        _mat[:,:,idx,2]=gen('{}/{}_Amp_zscore_2D_OK.csv'.format(manipPath,experiment),delimiter=',')

        #Get 1D Zscore
        zscore_1D_for_plot = np.nanmax(_mat[:,:,idx,2], axis=0)
            
    raw_zebrin_values = zebrins.loc['MEAN (normalized)'].values[:8]
    raw_zebrin_stds = zebrins.loc['STD (normalized)'].values[:8]
 
  
#FOR 2D ANALYSIS--------------------------------------------------------------------------------------------                
  
    fig, ax = plt.subplots(H,1, figsize=(7,9),sharex=True,sharey=True)
    plt.suptitle('{} 2D maps line by line'.format(condition))
    
    _MEDIAN_ZSCORE_2D, _AVERAGE_AMP_2D, _COUNT_2D, _POSITIONS_2D, _SUM_2D = [], [],[],[],[]
    
    for j in range(H):
        for y in range(N):

            #Create basis for concatenation at first loop
            if y == 0 :
                POSITIONS_2D = _mat[j,:,y,1]
                ZSCORES_2D = _mat[j,:,y,2]
                AMPS_2D = _mat[j,:,y,0]
            
            #Concatenate patterns for next loops
            else :
                POSITIONS_2D = np.concatenate((POSITIONS_2D,_mat[j,:,y,1]),axis=0)
                ZSCORES_2D = np.concatenate((ZSCORES_2D,_mat[j,:,y,2]),axis=0)
                AMPS_2D = np.concatenate((AMPS_2D,_mat[j,:,y,0]), axis=0)
                
            #SORT AMPLS AND ZSCORE ACCORDING TO POSITIONS
        SORTED_2D_AMPS = [x for _, x in sorted(zip(POSITIONS_2D,AMPS_2D))]
        SORTED_2D_ZSCORES = [x for _, x in sorted(zip(POSITIONS_2D,ZSCORES_2D))]
        
        ax[j].plot(sorted(POSITIONS_2D),SORTED_2D_ZSCORES,color=colors[i])
        ax[j].plot(sorted(POSITIONS_2D),np.ones(len(sorted(POSITIONS_2D)))*zscoreLimit,linestyle='--')
        label_line = j+1
        ax[j].set_ylabel('Zscore line {}'.format(label_line))
        
        if j == H-1:
            ax[j].set_xlabel('Distance (P1- norm)')
                
        zebrin_color='green'
        
        #P2+ contra zebrin band
        ax[j].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)        
        #P1+ zebrin band
        ax[j].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)       
        #P2+ ipsi zebrin band
        ax[j].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)     
        #PC_cluster
        ax[j].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymax=0.1, alpha=0.5,color='red')

        #BINNING FOR MEDIAN CALCUL
        step = binForMedian #In % of P1- : 10 is the last used 
        binning = np.arange(-left,right+step,step)
        
        _MEDS, _MADS, _POS, _COUNTS, _AMPS, _SUM = [],[],[],[],[],[]
        
        for y in range(len(binning)):
            
            if y == len(binning)-1:
                break
            
            start, stop = binning[y],binning[y+1]
            _meds, _mads, _pos, _count, _amps, _sum = [],[],[],[],[],[]
            
            #print ('Bin %s to %s'%(start, stop))
            
            SORTED_POSITIONS = sorted(POSITIONS_2D)
            
            for j in range(len(SORTED_POSITIONS)):
                if start < SORTED_POSITIONS[j] <= stop:
                    if np.isnan(SORTED_2D_ZSCORES[j])==False:
                        _meds.append(SORTED_2D_ZSCORES[j])
                        _pos.append(SORTED_POSITIONS[j])
                        _amps.append(SORTED_2D_AMPS[j])
                        _sum.append(SORTED_2D_AMPS[j])
                    
            _MEDS.append(np.nanmedian(_meds))
            _COUNTS.append(np.count_nonzero(_meds))
            _POS.append(np.nanmedian(_pos))
            _MADS.append(MAD(_meds, axis=0))
            _AMPS.append(np.nanmean(_amps,axis=0))
            _SUM.append(np.nansum(_sum,axis=0))
        
        _MEDIAN_ZSCORE_2D.append(np.asarray(_MEDS))
        _AVERAGE_AMP_2D.append(np.asarray(_AMPS))
        _COUNT_2D.append(np.asarray(_COUNTS))
        _POSITIONS_2D.append(np.asarray(_POS))
        _SUM_2D.append(np.asarray(_SUM))

    if saveFig == True:
        plt.savefig('{}/{}_2D_raw_sorting.pdf'.format(outputDir,condition)) 
        plt.savefig('{}/{}_2D_raw_sorting.png'.format(outputDir,condition))         


    fig, ax = plt.subplots(2,1,figsize=(14,5))

    plt.suptitle('{} 2D maps'.format(condition))
    ax[0].set_title('Median Zscore')
    median_zscore_2d = ax[0].imshow(stack_lines(_MEDIAN_ZSCORE_2D),interpolation=interpolationType, cmap='magma',vmin=minz,vmax=maxz,aspect='auto')
    fig.colorbar(median_zscore_2d, ax=ax[0])
        
    ax[1].set_title('Average Amplitude')
    ax[1].set_xticks(np.arange(0,len(_POSITIONS_2D[0]),1))
    ax[1].set_xticklabels(_POSITIONS_2D[0].astype(int),rotation=-90)
    mean_amplitude_2d = ax[1].imshow(stack_lines(_AVERAGE_AMP_2D),interpolation=interpolationType, cmap= 'magma_r',vmax=vmax,vmin=vmin,aspect='auto')
    fig.colorbar(mean_amplitude_2d,ax=ax[1])   

    if saveFig == True:    
        plt.savefig('{}/{}_2D_MedianZscore_and_AvgAmplitude.pdf'.format(outputDir,condition)) 
        plt.savefig('{}/{}_2D_MedianZscore_and_AvgAmplitude.png'.format(outputDir,condition))   
        
    if saveData == True:
        np.savetxt('{}/{}_2D_MedianZscore.csv'.format(outputDir,condition),
                   stack_lines(_MEDIAN_ZSCORE_2D),delimiter=',')
        
        np.savetxt('{}/{}_2D_AvgAmplitude.csv'.format(outputDir,condition),
                   stack_lines(_AVERAGE_AMP_2D),delimiter=',')

        np.savetxt('{}/{}_POSITIONAL_ARRAY.csv'.format(outputDir,condition),
                   binning,delimiter=',')      

        
    