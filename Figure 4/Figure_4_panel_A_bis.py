# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:33:53 2020

@author: ludov
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
import pingouin
import pandas as pd
import seaborn as sn
from statannot import add_stat_annotation
import sys
from scipy.stats import median_absolute_deviation as MAD

#Define the groups
groups = ['WT','ENR1','ENR2','LC','LS','EC','ES']
colors = ['skyblue','limegreen','green','lightcoral','black','orange','purple']
pal = ['skyblue','limegreen','green','lightcoral','grey','orange','purple']

#Pairs to analyse for stats
pairs= [("WT", "ENR1"), ("WT", "ENR2"),("WT", "LC"), ("WT", "LS"), ("WT", "EC"),("WT", "ES")]


#Input folder
dataSource = 'E:/000_PAPER/Amplitude_Analysis/Sigma'

#Savedir 
saveDir =  'E:/000_PAPER/Amplitude_Analysis/Sigma/03_SYNAPTIC/01_STATS'

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'


#Redirect console output to a txt file instead of console 
sys.stdout = open('{}/Activation_Statistics.txt'.format(saveDir),'w')


zscoreCut = 3.09

test = 'Mann-Whitney'
testB = 'Mann-Whitney'

#Apply stat correction ? 'bonferroni' or None
correction = None

bins = 200

#Do we save the data ?
saveData = True

#Do we save the figure ?
saveFig = False

#Do we perform the stats ?
statsToDo = True

fig, ax = plt.subplots(1,3, figsize=(10,9))
ax[0].set_ylabel('Amplitude (pA+/-SEM)')
ax[0].set_title('Avg synaptic Amplitude')
ax[1].set_ylabel('Total Amplitude (pA+/-SEM)')
ax[1].set_title('Total synaptic Amplitude')
ax[2].set_ylabel('Zscore > {} (%+/-SEM)'.format(zscoreCut))
ax[2].set_title('Proportion of active sites')
plt.tight_layout()

fig2,axx = plt.subplots(len(groups),1, figsize=(5,9), sharex=True, sharey=True)
plt.tight_layout()

fig3,axxx = plt.subplots(1,1)

fig4, box = plt.subplots(1,3)


for i in range(3):
    ax[i].set_xticks(range(len(groups)))
    ax[i].set_xticklabels(groups)
    
def SEM(data):
    return np.nanstd(data)/np.sqrt(len(data))

#FOR THE WHOLE MAP
avg, total, prop = [],[],[]

MEASURES = []

for group,index in zip(groups,range(len(groups))):
    print('Group = {}'.format(group))
    
    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir)]
        
    
    AVG_MEASURE, TOTAL_MEASURE, AVG_PROP = [],[],[]
    for manip,idx in zip(listOfExperiments,range(len(listOfExperiments))):
                
        #Get amplitudes or charges
        measures = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,fileType),delimiter=',').ravel())
        #Get corresponding zscores
        zscores = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,zscoreFileType),delimiter=',').ravel())

        #Extract significant measures on the whole map 
        sitesToKeep = [x for (x,y) in zip(measures,zscores) if y >= zscoreCut]
        
        #Remove NaNs from measures to compute proper proportion of active sites
        cleanedMeasures = [x for x in measures if np.isnan(x)==False]
        
        #Avg signal and proportion of active sites
        avgMeasure = round(np.nanmedian(sitesToKeep),3)
        proportionActiveSite = round(len(sitesToKeep)/len(cleanedMeasures)*100,3)
        totalMeasure = round(np.nansum(sitesToKeep),3)
        
        #Store distribution
        if idx == 0:
            distributions = sitesToKeep
        else:
            distributions= np.concatenate((distributions,sitesToKeep))
        
        #Append in list
        AVG_MEASURE.append(avgMeasure)
        TOTAL_MEASURE.append(totalMeasure)
        AVG_PROP.append(proportionActiveSite)
        
        
#        if group == 'LS':
#            box[0].annotate(manip,(index,avgMeasure))
#            box[1].annotate(manip,(index,totalMeasure))
#            box[2].annotate(manip,(index,proportionActiveSite))
            
     
    print ('Avg Amp+/-SD (map wise): {}+/-{}'.format(round(np.nanmean(AVG_MEASURE),2), round(np.nanstd(AVG_MEASURE),2)))
    print ('Avg activation (%)+/-SD: {}+/-{}'.format(round(np.nanmean(AVG_PROP),2), round(np.nanstd(AVG_PROP),2)))
    print ('Med activation (%)+/-SD: {}+/-{}'.format(round(np.nanmedian(AVG_PROP),2), round(MAD(AVG_PROP),2)))

            
    #Scatter plot             
    ax[0].scatter(np.ones(len(AVG_MEASURE))*index,AVG_MEASURE,color=colors[index])
    ax[0].bar(index,np.nanmean(AVG_MEASURE),alpha=0.5,color=colors[index])
    ax[0].errorbar(index,np.nanmean(AVG_MEASURE),yerr=SEM(AVG_MEASURE),color=colors[index])
    
    avg.append(AVG_MEASURE)

    ax[1].scatter(np.ones(len(TOTAL_MEASURE))*index,TOTAL_MEASURE,color=colors[index])
    ax[1].bar(index,np.nanmean(TOTAL_MEASURE),alpha=0.5,color=colors[index])        
    ax[1].errorbar(index,np.nanmean(TOTAL_MEASURE),yerr=SEM(TOTAL_MEASURE),color=colors[index])
    
    total.append(TOTAL_MEASURE)
    
    ax[2].scatter(np.ones(len(AVG_PROP))*index,AVG_PROP,color=colors[index])
    ax[2].bar(index,np.nanmean(AVG_PROP),alpha=0.5,color=colors[index])
    ax[2].errorbar(index,np.nanmean(AVG_PROP),yerr=SEM(AVG_PROP),color=colors[index])
    
    prop.append(AVG_PROP)
    
    axx[index].hist(distributions,bins=bins,density=True,color=colors[index],label=group)    
    axx[index].set_ylabel('Norm. Occurence')
    axx[index].set_xlabel('Amplitude (pA)')
    axx[index].legend(loc='best')
    
    axxx.hist(distributions,bins=bins,density=True,cumulative=True,color=colors[index],
              histtype='step',label=group)
    
    axxx.set_ylabel('Norm. count')
    axxx.set_xlabel('Amplitude (pA)')
    axxx.legend(loc='best')
    
    #Create a dataframe with all the distributions  
    
    indexList = ['Map#','Condition','Average amplitudes (pA)','Total amplitudes (pA)','Propotion (%)']

    groupDf = pd.DataFrame(np.vstack((listOfExperiments,
                                      [group for x in range(len(listOfExperiments))],
                                      AVG_MEASURE,
                                      TOTAL_MEASURE,
                                      AVG_PROP)),index=indexList).transpose()
    
    #Pass numeric columns to float
    groupDf['Average amplitudes (pA)'] = groupDf['Average amplitudes (pA)'].astype(float)
    groupDf['Total amplitudes (pA)'] = groupDf['Total amplitudes (pA)'].astype(float)
    groupDf['Propotion (%)'] = groupDf['Propotion (%)'].astype(float)
    
    MEASURES.append(groupDf)
    



if saveFig == True:
    fig.savefig('{}/Amplitudes.pdf'.format(saveDir))
    fig.savefig('{}/Amplitudes.png'.format(saveDir))
    
    fig2.savefig('{}/Amplitudes_Distributions.pdf'.format(saveDir))
    fig2.savefig('{}/Amplitudes_Distributions.png'.format(saveDir))
    
    fig3.savefig('{}/Cumulative_Amplitudes_Distributions.pdf'.format(saveDir))
    fig3.savefig('{}/Cumulative_Amplitudes_Distributions.png'.format(saveDir))
    


#Create global dataframe for stats and shit
allData = pd.concat(MEASURES)


#Fill the box plot 
sn.boxplot(x=allData["Condition"],y=allData["Average amplitudes (pA)"],ax=box[0],palette=pal)
sn.swarmplot(x=allData["Condition"],y=allData["Average amplitudes (pA)"],ax=box[0],color='0.25')

add_stat_annotation(box[0], x=allData["Condition"], y=allData["Average amplitudes (pA)"], order=groups,
                    box_pairs=pairs,
                    test=test, text_format='star', loc='inside', verbose=2,comparisons_correction=correction)


sn.boxplot(x=allData["Condition"],y=allData["Total amplitudes (pA)"],ax=box[1],palette=pal)
sn.swarmplot(x=allData["Condition"],y=allData["Total amplitudes (pA)"],ax=box[1],color='0.25')

add_stat_annotation(box[1], x=allData["Condition"], y=allData["Total amplitudes (pA)"], order=groups,
                    box_pairs=pairs,
                    test=test, text_format='star', loc='inside', verbose=2,comparisons_correction=correction)

sn.boxplot(x=allData["Condition"],y=allData["Propotion (%)"],ax=box[2],palette=pal)
sn.swarmplot(x=allData["Condition"],y=allData["Propotion (%)"],ax=box[2],color='0.25')

add_stat_annotation(box[2], x=allData["Condition"], y=allData["Propotion (%)"], order=groups,
                    box_pairs=pairs,
                    test=testB, text_format='star', loc='inside', verbose=2,comparisons_correction=correction)

if saveData == True:
    allData.to_excel('{}/Average_Amplitudes.xlsx'.format(saveDir))

if statsToDo == True:

    kruskalAvg = pingouin.kruskal(data = allData, dv='Average amplitudes (pA)',between='Condition')
    print ('Avg Amp')
    print (kruskalAvg)
    
    kruskalTotal = pingouin.kruskal(data = allData, dv='Total amplitudes (pA)',between='Condition')
    print ('Total Amp')
    print (kruskalTotal)
    
    kruskalProportion = pingouin.kruskal(data = allData, dv='Propotion (%)',between='Condition')
    print ('Propotion of active sites')
    print (kruskalProportion)
    
    
    anovaAvg = pingouin.anova(data = allData, dv='Average amplitudes (pA)',between='Condition')
    print ('Avg Amp')
    print (anovaAvg)
    
    anovaTotal = pingouin.anova(data = allData, dv='Total amplitudes (pA)',between='Condition')
    print ('Total Amp')
    print (anovaTotal)
    
    anovaProportion = pingouin.anova(data = allData, dv='Propotion (%)',between='Condition')
    print ('Propotion of active sites')
    print (anovaProportion)
