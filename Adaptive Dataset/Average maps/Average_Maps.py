# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:22:58 2020

This code computes median synaptic patterns and median maps from individual synaptic maps, as well as 2D maps and zscores
AMPLITUDE way : inputs are amplitude (pA) based maps 

@author: ludov
"""

import numpy as np 
from numpy import genfromtxt as gen
from matplotlib import pyplot as plt 
import os 
import pandas as pd 
import matplotlib
from scipy import stats
import seaborn as sn 
matplotlib.rcParams['pdf.fonttype'] = 42

#The groups to analyse
conditions = ['WT','ENR1','ENR2','EC','ES','LC','LS']
colors = ['skyblue','limegreen','green','orange','grey','lightcoral','black']

#The sheetnames in excel file
sheets=['WT','ENR','ENR','EC','ES','LC','LS']


#General directory to find the data
inputDir = 'D:/03_FORMATED_DATA/For Paper/EPHYS/Adaptive_Dataset'
#inputDir = './For Paper/EPHYS/Adaptive_Dataset'

#Where to save datas/figures
#outputPath = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/LTP-LTD question/code output'
outputPath = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/LTP-LTD question/code output'
#outputPath = './data/'
#Zebrin file
zebrinFile = 'D:/03_FORMATED_DATA/For Paper/EPHYS/Mesures_ZII_HighRes_WT_ENR_EC_LC_ES_LS.xlsx'
#zebrinFile = './For Paper/EPHYS/Mesures_ZII_HighRes_WT_ENR_EC_LC_ES_LS.xlsx'

#To constrict 1D plots
ylim = 20

binForMedian = 10  #In % of P1- : 10 is the last value used

left,right = 210,210 #Borders in %P1- 

zscoreLimit =2.0 #Limit of significance for zscore

SDfactor = 2

vmax = 75
minz, maxz = 2,4

#Do we save anything ?
saveFig = False #For figures
saveData = False

varFig, varax = plt.subplots(1,1)

#------------------------------------FUNCTIONS---------------------------------------------------------
#---------------------------------DO NOT MODIFY--------------------------------------------------------
def stack_lines(_list):
    
    stacked_map = np.vstack((_list[0],
                            _list[1],
                            _list[2],
                            _list[3],
                            _list[4],
                            _list[5]))
    
    return stacked_map


averagedMaps = []

#Iterate for each condition
for condition,sheet,i in zip(conditions,sheets,range(len(conditions))): 
    
    parentDir = '{}/{}'.format(inputDir,condition)
    print(parentDir)
    
    listOfExperiments = os.listdir(parentDir)
    
    #Load zebrin file (i.e. corresponding datasheet)
    zebrins = pd.read_excel(zebrinFile,sheet_name=sheet,
                            index_col=0,header=1)
    
    #Matrix to append maps and positions
    H = 6           #Map height in sites
    L = 64          #Map width in sites
    N = len(listOfExperiments)
    
    _mat = np.zeros((H,L,N,3)) #[0] for map, [1] for position and [2] for Zscore

    #Get the noise values to establish threshold for 2D maps
    noiseValues=[]

    #Iterate on each experiment
    for experiment,idx in zip(listOfExperiments,range(N)):
        #print (experiment)        
        manipPath = '{}/{}'.format(parentDir,experiment)
        
        #Load map in Matrix
        _mat[:,:,idx,0]=gen('{}/{}_Amp_2D_OK.csv'.format(manipPath,experiment),delimiter=',')
        
        #Get the positions
        pos = gen('{}/{}_Positions_cp_centered_OK.csv'.format(manipPath,experiment),delimiter=',')                
        pos_2D = (pos,pos,pos,pos,pos,pos)        
        _mat[:,:,idx,1]=np.reshape(pos_2D,(H,L))
        
        #And now the 2D Zscore    
        _mat[:,:,idx,2]=gen('{}/{}_Amp_zscore_2D_OK.csv'.format(manipPath,experiment),delimiter=',')
        
        noiseValues.append(np.nanmean(np.abs(gen('{}/{}_Amp_Noisemap_OK.csv'.format(manipPath,experiment),delimiter=',').ravel())))
        
#FOR 2D ANALYSIS--------------------------------------------------------------------------------------------
#FOR 2D ANALYSIS--------------same shit but line by line first----------------------------------------------        
#FOR 2D ANALYSIS--------------------------------------------------------------------------------------------                
  
    _MEDIAN_ZSCORE_2D, _AVERAGE_AMP_2D, _COUNT_2D, _POSITIONS_2D, _SUM_2D, _MEDIAN_2D, _VAR2D = [],[],[],[],[],[],[]
    
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
        
        label_line = j+1
        
                
        zebrin_color='green'
        
        raw_zebrin_values = zebrins.loc['MEAN (normalized)'].values[:8]
        raw_zebrin_stds = zebrins.loc['STD (normalized)'].values[:8]
        
   
        #BINNING FOR MEDIAN CALCUL
        step = binForMedian #In % of P1- : 10 is the last used 
        binning = np.arange(-left,right+step,step)
        
        _MEDS, _MADS, _POS, _COUNTS, _AMPS, _SUM, _AMPMED, _VAR= [],[],[],[],[],[],[],[]
        
        for y in range(len(binning)):
            
            if y == len(binning)-1:
                break
            
            start, stop = binning[y],binning[y+1]
            _meds, _mads, _pos, _count, _amps, _sum, _ampmed = [],[],[],[],[],[],[]
            
            #print ('Bin %s to %s'%(start, stop))
            
            SORTED_POSITIONS = sorted(POSITIONS_2D)
            
            for j in range(len(SORTED_POSITIONS)):
                if start < SORTED_POSITIONS[j] <= stop:
                    if np.isnan(SORTED_2D_ZSCORES[j])==False:
                        _meds.append(SORTED_2D_ZSCORES[j])
                        _pos.append(SORTED_POSITIONS[j])
                        _amps.append(SORTED_2D_AMPS[j])
                        _sum.append(SORTED_2D_AMPS[j])
                        _ampmed.append(SORTED_2D_AMPS[j])
                    
            _MEDS.append(np.nanmedian(_meds))
            _COUNTS.append(np.count_nonzero(_meds))
            _POS.append(np.nanmedian(_pos))
            _AMPS.append(np.nanmean(_amps,axis=0))
            _SUM.append(np.nansum(_sum,axis=0))
            _AMPMED.append(np.nanmedian(_ampmed))
            _VAR.append(stats.variation(_amps))
        
        _MEDIAN_ZSCORE_2D.append(np.asarray(_MEDS))
        _AVERAGE_AMP_2D.append(np.asarray(_AMPS))
        _COUNT_2D.append(np.asarray(_COUNTS))
        _POSITIONS_2D.append(np.asarray(_POS))
        _SUM_2D.append(np.asarray(_SUM))
        _MEDIAN_2D.append(np.asarray(_AMPMED))
        _VAR2D.append(np.asarray(_VAR))
  
    fig, ax = plt.subplots(3,1,figsize=(14,8))

    vmin = np.nanmean(noiseValues) + SDfactor * np.nanstd(noiseValues)

    plt.suptitle('{} 2D maps'.format(condition))
    ax[0].set_title('Median Zscore')
    median_zscore_2d = ax[0].imshow(stack_lines(_MEDIAN_ZSCORE_2D),interpolation='sinc', cmap='magma',vmin=minz,vmax=maxz,aspect='auto')
    fig.colorbar(median_zscore_2d, ax=ax[0])
        
    ax[1].set_title('Average Amplitude')
    # ax[1].set_xticks(np.arange(0,len(_POSITIONS_2D[0]),1))
    # ax[1].set_xticklabels(_POSITIONS_2D[0].astype(int),rotation=-90)
    mean_amplitude_2d = ax[1].imshow(np.abs(stack_lines(_AVERAGE_AMP_2D)),interpolation='sinc', cmap= 'magma',vmax=vmax,vmin=vmin,aspect='auto')
    fig.colorbar(mean_amplitude_2d,ax=ax[1])   
    
    averagedMaps.append(np.abs(stack_lines(_AVERAGE_AMP_2D)))
    
    ax[2].set_title('AMP CV')
    ax[2].set_xticks(np.arange(0,len(_POSITIONS_2D[0]),1))
    ax[2].set_xticklabels(_POSITIONS_2D[0].astype(int),rotation=-90)
    fano = ax[2].imshow(np.abs(stack_lines(_VAR2D)),interpolation='sinc', cmap= 'magma',aspect='auto',vmax=2)
    fig.colorbar(fano, ax=ax[2])
    
    sn.kdeplot(np.abs(_VAR2D).ravel(), ax=varax,label=condition,color=colors[i])
    varax.legend(loc='best')
    varax.set_title('Amp CV density distribution')
    
    
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
        
        np.savetxt('{}/{}_2D_MedianAmplitude.csv'.format(outputPath,condition),
                   stack_lines(_MEDIAN_2D),delimiter=',')
        
        
#Compare 2D maps 
compFig, compAx = plt.subplots(len(conditions),1,figsize=(5,9))
compChart, compPie = plt.subplots(len(conditions),1,figsize=(5,9))

for i in range(len(averagedMaps)):
    
    #Divide maps by control
    deltaMap = averagedMaps[i] / averagedMaps[0]
    
    ltp = len([x for x in np.nditer(deltaMap) if x>1])/len(deltaMap.ravel())*100
    ltd = len([x for x in np.nditer(deltaMap) if x<1])/len(deltaMap.ravel())*100
    
    compPie[i].pie([ltp,ltd],colors=['crimson','royalblue'], labels=['LTP','LTD'],autopct='%.0f%%')

    deltaPlot = compAx[i].imshow(deltaMap,cmap='coolwarm',interpolation='sinc', vmin=0, vmax=2)
    compAx[i].set_title('(WT/{})'.format(conditions[i]))
    compFig.colorbar(deltaPlot,ax=compAx[i])   
    
    
    pd.DataFrame(data=averagedMaps[i]).to_excel('{}/{}_Average2Dmap.xlsx'.format(outputPath,conditions[i]))
    
    

