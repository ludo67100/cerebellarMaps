# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 12:37:26 2020

@author: Ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
from scipy import stats 



#Define the groups
groups = ['P9P10','P12P13','P14P18','P30P40']
colors = ['lightskyblue','skyblue','deepskyblue','royalblue']
pal = ['lightskyblue','skyblue','deepskyblue','royalblue']

#Pairs to analyse for stats
pairs= [("P9P10", "P12P13"), ("P9P10", "P14P18"), ("P9P10", "P30P40"), ("P12P13", "P14P18"),
        ("P12P13", "P30P40"),("P14P18", "P30P40")]


#Input folder
dataSource = 'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS'

#Savedir 
saveDir =  'E:/000_PAPER/Development/Amplitude_Analysis/03_SYNAPTIC'

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'



zscoreCut = 3.09

test = 'Mann-Whitney'
testB = 'Mann-Whitney'

#Apply stat correction ? 'bonferroni' or None
correction = 'bonferroni'

bins = 20

loclimitPos = 33.00
loclimitNeg = -0.0

#Do we save the data ?
saveData = True

#Do we save the figure ?
saveFig = True

#Do we perform the stats ?
statsToDo = True

fig, ax = plt.subplots(1,2)
ax[0].set_ylabel('Avg Amplitude (pA) +/-SEM')
ax[1].set_ylabel('% of active sites +/-SEM')
    
def SEM(data):
    return np.nanstd(data)/np.sqrt(len(data))

#FOR THE WHOLE MAP

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
        
        ax[0].plot([index-0.25, index+0.25], [individualLocal[i],individualDist[i]], color='0.5', marker='x',alpha=0.3)
        
        ax[1].plot([index-0.25, index+0.25], [LocalCount[i],DistCount[i]], color='0.5', marker='x',alpha=0.3)
            
            
    avgLocal = np.nanmean(individualLocal)
    stdLocal = np.nanstd(individualLocal)
    semLocal = SEM(individualLocal)
    
    avgDist = np.nanmean(individualDist)
    stdDist = np.nanstd(individualDist)
    semDist = SEM(individualDist)
    
    
    countLocalavg = np.nanmean(LocalCount)
    countLocalstd = np.nanstd(LocalCount)
    countLocalSEM = SEM(LocalCount)
    
    countDistavg = np.nanmean(DistCount)
    countDiststd = np.nanstd(DistCount)
    countDistSEM = SEM(DistCount)
    
    ax[0].plot([index-0.25, index+0.25], [avgLocal,avgDist], color=colors[index], marker='o',label=group)
    ax[0].fill_between([index-0.25, index+0.25], [avgLocal-semLocal,avgDist-semDist],[avgLocal+semLocal,avgDist+semDist],color=colors[index],alpha=0.3)
        
    ax[1].plot([index-0.25, index+0.25], [countLocalavg,countDistavg], color=colors[index], marker='o',label=group)
    ax[1].fill_between([index-0.25, index+0.25], [countLocalavg-countLocalSEM,countDistavg-countDistSEM],
                       [countLocalavg+countLocalSEM,countDistavg+countDistSEM],color=colors[index],alpha=0.3)
    
    
    
    #Wilcoxon test for amplitudes 
    statAmp , pAmp = stats.wilcoxon(individualLocal,individualDist)
    statProp, pProp = stats.wilcoxon(LocalCount,DistCount)
    
    print ('Wilcoxon signed rank test')
    print ('Loc vs Dist amplitudes p = {}'.format(pAmp))
    print ('Lov vs Dist Proportion p = {}'.format(pProp))
    

ax[0].legend(loc='best')
ax[1].legend(loc='best')
    
    
        
            
                
                
        
        
        
        