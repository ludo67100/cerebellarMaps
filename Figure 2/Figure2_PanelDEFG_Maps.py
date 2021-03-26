# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:22:58 2020

@author: ludov
"""

#--------------FILL HERE INPUT/OUTPUT-------------------------------------------
#Input folder : please put path up to EPHYS dataset (e.g., 'C/USER/FOLDER/EPHYS')
inputDir = 'E:/03_FORMATED_DATA/For Paper/EPHYS'

#Savedir  : please specify a path where to save the output
saveDir =  'E:/code_output_test'
#-------------------------------------------------------------------------------

dataSource = '{}/Development_Dataset'.format(inputDir)

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

#Zebrin file
zebrinFile = '{}/Mesures_ZII_LowRes_Adult_and_Dev.xlsx'.format(inputDir)

#To constrict 1D plots
ylim = 20

binForMedian = 10  #In % of P1- : 10 is the last value used

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
    
    parentDir = '{}/{}'.format(dataSource,condition)
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

    #First plot: all patterns
    all_maps, maps = plt.subplots(N,1,sharex=True,sharey=True,figsize=(8,9))
    all_maps.suptitle('{} Zscore (amplitude) patterns'.format(condition))
    

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
        
        #Plot the zscore patterns
        maps[idx].fill_between(pos,0,zscore_1D_for_plot,label='{}'.format(experiment),color=colors[i])
        maps[idx].legend(loc='best')

        #FOR 1D ANALYSIS-------------------------------------------------------
        #Create basis for concatenation at first loop
        if idx == 0 :
            POSITIONS_1D = pos
            ZSCORES_1D = np.max(_mat[:,:,idx,2], axis=0)
            AMPS_1D = np.min(_mat[:,:,idx,0], axis=0)
        
        #Concatenate patterns for next loops
        else :
            POSITIONS_1D = np.concatenate((POSITIONS_1D,pos),axis=0)
            ZSCORES_1D = np.concatenate((ZSCORES_1D,np.max(_mat[:,:,idx,2], axis=0)),axis=0)
            AMPS_1D = np.concatenate((AMPS_1D,np.min(_mat[:,:,idx,0], axis=0)), axis=0)
            
    #Save first figure
    if saveFig == True:
        #Create output directory
        outputPath = '{}/{}'.format(outputDir,condition)
        try:
            os.makedirs(outputPath)
        except:
            pass
        
        plt.savefig('{}/{}_all_1D_zscores.pdf'.format(outputPath,condition))
        plt.savefig('{}/{}_all_1D_zscores.png'.format(outputPath,condition))

        
    #SORT AMPLS AND ZSCORE ACCORDING TO POSITIONS
    SORTED_1D_AMPS = [x for _, x in sorted(zip(POSITIONS_1D,AMPS_1D))]
    SORTED_1D_ZSCORES = [x for _, x in sorted(zip(POSITIONS_1D,ZSCORES_1D))]
            
    #Then plot
    fig, ax = plt.subplots(2,1, figsize=(7,8))
    plt.suptitle('{} raw sorting (max values per maps)'.format(condition))
    
    ax[0].plot(sorted(POSITIONS_1D),SORTED_1D_AMPS, label='Amp',color=colors[i])
    ax[0].set_ylabel('Amp (pA)')
    
    ax[1].plot(sorted(POSITIONS_1D),np.ones(len(sorted(POSITIONS_1D)))*2.0, color='gray',linestyle='--')
    ax[1].plot(sorted(POSITIONS_1D),np.ones(len(sorted(POSITIONS_1D)))*3.09, color='black',linestyle='--')

    ax[1].plot(sorted(POSITIONS_1D),SORTED_1D_ZSCORES, label='Max Zscores',color=colors[i])
    ax[1].set_ylabel('Zscore')
    ax[1].set_xlabel('Position (P1- norm)')
    
    raw_zebrin_values = zebrins.loc['MEAN (normalized)'].values[:8]
    raw_zebrin_stds = zebrins.loc['STD (normalized)'].values[:8]

    zebrin_color='green'
    
    #P2+ contra zebrin band
    ax[0].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)
    ax[1].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)
    
    #P1+ zebrin band
    ax[0].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)
    ax[1].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)
    
    #P2+ ipsi zebrin band
    ax[0].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)
    ax[1].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)
    
    #PC_cluster
    ax[0].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymin=0.9, alpha=0.5,color='red')
    ax[1].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymax=0.1, alpha=0.5,color='red')
   
    if saveData == True: 
        #Create output directory
        outputPath = '{}/{}'.format(saveDir,condition)
        try:
            os.makedirs(outputPath)
        except:
            pass
        
        np.savetxt('{}/{}_SORTED_1D_AMPS.csv'.format(outputPath,condition),np.asarray(SORTED_1D_AMPS),delimiter=',')
        np.savetxt('{}/{}_SORTED_1D_ZSCORE.csv'.format(outputPath,condition),np.asarray(SORTED_1D_ZSCORES),delimiter=',')
        np.savetxt('{}/{}_SORTED_1D_POSITIONS.csv'.format(outputPath,condition),np.asarray(sorted(POSITIONS_1D)),delimiter=',')


    if saveFig == True:        
        plt.savefig('{}/{}_1D_raw_sorting.pdf'.format(outputPath,condition))
        plt.savefig('{}/{}_1D_raw_sorting.png'.format(outputPath,condition))
                
    #BINNING FOR MEDIAN CALCUL
    step = binForMedian #In % of P1- : 12 c'est pas mal 
    binning = np.arange(-left,right+step,step)
    
    _MEDS, _MADS, _POS, _COUNTS = [],[],[],[]
    
    for y in range(len(binning)):
        
        if y == len(binning)-1:
            break
        
        start, stop = binning[y],binning[y+1]
        _meds, _mads, _pos, _count = [],[],[],[]
                
        SORTED_POSITIONS = sorted(POSITIONS_1D)
        
        for j in range(len(SORTED_POSITIONS)):
            if start < SORTED_POSITIONS[j] <= stop:
                if np.isnan(SORTED_1D_ZSCORES[j])==False:
                    _meds.append(SORTED_1D_ZSCORES[j])
                    _pos.append(SORTED_POSITIONS[j])
                
        _MEDS.append(np.nanmedian(_meds))
        _COUNTS.append(np.count_nonzero(_meds))
        _POS.append(np.nanmedian(_pos))
        _MADS.append(MAD(_meds, axis=0))
        
    #Save the data
    if saveData == True: 
        np.savetxt('{}/{}_MEDIAN_ZSCORE_1D_{}_BINS.csv'.format(outputPath,condition,binForMedian),np.asarray(_MEDS),delimiter=',')
        np.savetxt('{}/{}_MEDIAN_MADS_1D_{}_BINS.csv'.format(outputPath,condition,binForMedian),np.asarray(_MADS),delimiter=',')
        np.savetxt('{}/{}_POSITIONS_1D_{}_BINS.csv'.format(outputPath,condition,binForMedian),np.asarray(_POS),delimiter=',')
                

#FOR 2D ANALYSIS--------------------------------------------------------------------------------------------
#FOR 2D ANALYSIS--------------same shit but line by line first----------------------------------------------        
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
        plt.savefig('{}/{}_2D_raw_sorting.pdf'.format(outputPath,condition)) 
        plt.savefig('{}/{}_2D_raw_sorting.png'.format(outputPath,condition))         


    fig, ax = plt.subplots(2,1,figsize=(14,5))

    plt.suptitle('{} 2D maps'.format(condition))
    ax[0].set_title('Median Zscore')
    median_zscore_2d = ax[0].imshow(stack_lines(_MEDIAN_ZSCORE_2D),interpolation=interpolationType, cmap='hot',vmin=minz,vmax=maxz,aspect='auto')
    fig.colorbar(median_zscore_2d, ax=ax[0])
        
    ax[1].set_title('Average Amplitude')
    ax[1].set_xticks(np.arange(0,len(_POSITIONS_2D[0]),1))
    ax[1].set_xticklabels(_POSITIONS_2D[0].astype(int),rotation=-90)
    mean_amplitude_2d = ax[1].imshow(stack_lines(_AVERAGE_AMP_2D),interpolation=interpolationType, cmap= 'hot_r',vmax=vmax,vmin=vmin,aspect='auto')
    fig.colorbar(mean_amplitude_2d,ax=ax[1])   

    if saveFig == True:    
        plt.savefig('{}/{}_2D_MedianZscore_and_AvgAmplitude.pdf'.format(outputPath,condition)) 
        plt.savefig('{}/{}_2D_MedianZscore_and_AvgAmplitude.png'.format(outputPath,condition))   
        
    if saveData == True:
        np.savetxt('{}/{}_2D_MedianZscore.csv'.format(outputPath,condition),
                   stack_lines(_MEDIAN_ZSCORE_2D),delimiter=',')
        
        np.savetxt('{}/{}_2D_AvgAmplitude.csv'.format(outputPath,condition),
                   stack_lines(_AVERAGE_AMP_2D),delimiter=',')

        np.savetxt('{}/{}_POSITIONAL_ARRAY.csv'.format(outputPath,condition),
                   binning,delimiter=',')   