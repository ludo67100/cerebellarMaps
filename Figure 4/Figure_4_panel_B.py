# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:55:49 2020

@author: Ludovic.spaeth
"""


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
from scipy import stats 
import pandas as pd
import seaborn as sn
import sys
from scipy.stats import median_absolute_deviation as MAD
import math 

#Define the groups
groups = ['WT','ENR1','ENR2','LC','LS','EC','ES']
colors = ['skyblue','limegreen','green','lightcoral','black','orange','purple']
pal = ['skyblue','limegreen','green','lightcoral','grey','orange','purple']


#General directory to find the data
dataSource = 'E:/000_PAPER/Amplitude_Analysis/Sigma'

saveDir = 'E:/000_PAPER/Amplitude_Analysis/Sigma/03_SYNAPTIC/01_STATS'

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'

zscoreCut = 3.09

test = 'Mann-Whitney'
testB = 'Mann-Whitney'

#Apply stat correction ? 'bonferroni' or None
correction = 'bonferroni'

bins = 20

loclimitPos = 30.00
loclimitNeg = -10.00

#Do we save the data ?
saveData = True

#Do we save the figure ?
saveFig = False

#Do we perform the stats ?
statsToDo = True

#Round factor (decimal)
RF =2

if saveData == True:
    #Redirect console output to a txt file instead of console 
    sys.stdout = open('{}/Local_vs_Distal_Statistics.txt'.format(saveDir),'w')

AmpFig, AmpAx = plt.subplots(1,(len(groups)),sharex=True,sharey=True,figsize=(16,9))
AmpAx[0].set_ylabel('Avg Amplitude (pA)')

PropFig, PropAx = plt.subplots(1,(len(groups)),sharex=True,sharey=True,figsize=(16,9))
PropAx[0].set_ylabel('Active sites (%)')
    
def SEM(data):
    return np.nanstd(data)/np.sqrt(len(data))

#FOR THE WHOLE MAP
    
LOCAL_AMPS = []
DISTAL_AMPS = []
LOCAL_PROP = []
DISTAL_PROP = []


for group,index in zip(groups,range(len(groups))):
    print('Group = {}'.format(group))
    
    globLocal, globDist = [],[]
    
    individualLocal, individualDist = [],[]
    
    LocalCount, DistCount = [],[]
    
    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir)]
        
    
    AVG_MEASURE, TOTAL_MEASURE, AVG_PROP, VAR_INDEX = [],[],[],[]
    for manip,idx in zip(listOfExperiments,range(len(listOfExperiments))):
        
        tempLoc, tempDist = [],[]
        LocAmount, DistAmount = [],[]
                
        #Get amplitudes or charges
        measures = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,fileType),delimiter=','))
        #Get corresponding zscores
        zscores = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,zscoreFileType),delimiter=','))
        
        #Get positions 
        positions = np.genfromtxt('{}/{}/{}_Positions_cp_centered_OK.csv'.format(inputDir,manip,manip,zscoreFileType),delimiter=',')
        
        
        for j in range(len(positions)): 
                            
            if loclimitNeg < positions[j] <= loclimitPos:
                
                tempLoc.append([x for (x,y) in zip(measures[:,j], zscores[:,j]) if y >= zscoreCut])
                LocAmount.append([x for (x,y) in zip(measures[:,j], zscores[:,j])])
                
            else:
                
                tempDist.append([x for (x,y) in zip(measures[:,j], zscores[:,j]) if y >= zscoreCut])
                DistAmount.append([x for (x,y) in zip(measures[:,j], zscores[:,j])])
                
            
        globLocal.append([x for xs in tempLoc for x in xs])
        globDist.append([x for xs in tempDist for x in xs])
        
        individualLocal.append(np.nanmean([x for xs in tempLoc for x in xs]))
        individualDist.append(np.nanmean([x for xs in tempDist for x in xs]))
        
        LocalCount.append(len([x for xs in tempLoc for x in xs])/len([x for xs in LocAmount for x in xs])*100)
        DistCount.append(len([x for xs in tempDist for x in xs])/len([x for xs in DistAmount for x in xs])*100)
        
    for i in range(len(individualLocal)): 
        
        AmpAx[index].plot([0,1], [individualLocal[i],individualDist[i]], color='0.5', marker='x',alpha=0.2)
        
        PropAx[index].plot([0,1], [LocalCount[i],DistCount[i]], color='0.5', marker='x',alpha=0.2)
        
    
    #Create dataframe to hold data
    AmpDf = pd.DataFrame()
    AmpDf['Local'] = individualLocal
    AmpDf['Distal'] = individualDist
    
    PropDf = pd.DataFrame()
    PropDf['Local'] = LocalCount
    PropDf['Distal'] = DistCount
    
    #Append in global dataframes
    LOCAL_AMPS.append(pd.DataFrame(individualLocal))
    DISTAL_AMPS.append(pd.DataFrame(individualDist))
    LOCAL_PROP.append(pd.DataFrame(LocalCount))
    DISTAL_PROP.append(pd.DataFrame(DistCount))
    
            
    avgLocal = np.nanmean(individualLocal)
    stdLocal = np.nanstd(individualLocal)
    semLocal = SEM(individualLocal)
    
    medLocal = np.nanmedian(individualLocal)
    madLocal = MAD([x for x in individualLocal if math.isnan(x)==False])
    
    avgDist = np.nanmean(individualDist)
    stdDist = np.nanstd(individualDist)
    semDist = SEM(individualDist)
    
    medDist = np.nanmedian(individualDist)
    madDist = MAD([x for x in individualDist if math.isnan(x)==False])
    
    
    countLocalavg = np.nanmean(LocalCount)
    countLocalstd = np.nanstd(LocalCount)
    countLocalSEM = SEM(LocalCount)
    
    countLocalMed = np.nanmedian(LocalCount)
    countLocalMad = MAD([x for x in LocalCount if math.isnan(x)==False])
    
    countDistavg = np.nanmean(DistCount)
    countDiststd = np.nanstd(DistCount)
    countDistSEM = SEM(DistCount)
    
    countDistMed = np.nanmedian(DistCount)
    countDistMad = MAD([x for x in DistCount if math.isnan(x)==False])
    
    print('SYNAPTIC AMPLITUDE')
    print('Avg amp(pA) in Loc = {} ({})'.format(round(avgLocal,RF),round(stdLocal,RF)))
    print('Avg amp(pA) in Dist = {} ({})'.format(round(avgDist,RF),round(stdDist,RF)))
    
    print('Med amp(pA) in Loc = {} ({})'.format(round(medLocal,RF),round(madLocal,RF)))
    print('Med amp(pA) in Dist = {} ({})'.format(round(medDist,RF),round(madDist,RF)))
   
    print('ACTIVE SITES')
    print('Avg % in Loc = {} ({})'.format(round(countLocalavg,RF),round(countLocalstd,RF)))
    print('Avg % in Dist = {} ({})'.format(round(countDistavg,RF),round(countDiststd,RF)))    
    
    print('MEd % in Loc = {} ({})'.format(round(countLocalMed,RF),round(countLocalMad,RF)))
    print('MEd % in Dist = {} ({})'.format(round(countDistMed,RF),round(countDistMad,RF)))  
    
#    ax[0].plot([index-0.25, index+0.25], [avgLocal,avgDist], color=colors[index], marker='o',label=group)
#    ax[0].fill_between([index-0.25, index+0.25], [avgLocal-semLocal,avgDist-semDist],[avgLocal+semLocal,avgDist+semDist],color=colors[index],alpha=0.3)
#        
#    ax[1].plot([index-0.25, index+0.25], [countLocalavg,countDistavg], color=colors[index], marker='o',label=group)
#    ax[1].fill_between([index-0.25, index+0.25], [countLocalavg-countLocalSEM,countDistavg-countDistSEM],
#                       [countLocalavg+countLocalSEM,countDistavg+countDistSEM],color=colors[index],alpha=0.3)

    sn.boxplot(data=AmpDf, color=colors[index],ax=AmpAx[index])
    sn.boxplot(data=PropDf, color=colors[index],ax=PropAx[index])  
    
    
    #Wilcoxon test for amplitudes 
    statAmp , pAmp = stats.wilcoxon(individualLocal,individualDist)
    statProp, pProp = stats.wilcoxon(LocalCount,DistCount)
    
    print ('Wilcoxon signed rank test')
    print ('Loc vs Dist amplitudes w = {}'.format(statAmp))
    print ('Loc vs Dist amplitudes p = {}'.format(pAmp))
    print ('Lov vs Dist Proportion w = {}'.format(statProp))
    print ('Lov vs Dist Proportion p = {}'.format(pProp))
    
    print ('')
    
    
    AmpAx[index].set_title(group)
    PropAx[index].set_title(group)
    
#Compare local and distal from each group 
    
LOCAL_AMPS = pd.concat(LOCAL_AMPS,ignore_index=True,axis=1)
DISTAL_AMPS = pd.concat(DISTAL_AMPS,ignore_index=True,axis=1)
LOCAL_PROP = pd.concat(LOCAL_PROP,ignore_index=True,axis=1)
DISTAL_PROP = pd.concat(DISTAL_PROP,ignore_index=True,axis=1)

dfs = [LOCAL_AMPS, DISTAL_AMPS, LOCAL_PROP, DISTAL_PROP]
titles = ['Local average amplitude','Distal average amplitude','Local % of active sites','Distal % of active sites']
labels = ['Amplitude (pA)','Amplitude (pA)','Active sites (%)','Active sites (%)']

for table, title, label in zip(dfs, titles, labels): 
    table.columns=groups
    
    plt.figure(figsize=(6,4))
    sn.boxplot(data=table, palette=pal)
    plt.ylabel(label)
    plt.title(title)



    
    