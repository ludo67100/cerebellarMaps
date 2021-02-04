# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:01:23 2020

@author: Ludovic.spaeth
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
from scipy import stats 
import sys



#Define the groups
groups = ['P9P10','P12P13','P14P18','P30P40']
colors = ['lightskyblue','skyblue','deepskyblue','royalblue']
pal = ['lightskyblue','skyblue','deepskyblue','royalblue']

colors = ['black','orange','lightcoral','skyblue']
pal = ['black','orange','lightcoral','skyblue']

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

#Redirect console output to a txt file instead of console 
sys.stdout = open('{}/Statistics.txt'.format(saveDir),'w')

##Input folder
#dataSource = 'E:/000_PAPER/Charge_Analysis/Sigma'
#
##Savedir 
#saveDir =  'D:/000_PAPER/Charge_Analysis/Sigma/03_SYNAPTIC'
#
##Target file type
#fileType = 'Charge_2D_OK.csv'
#zscoreFileType = 'Charge_zscore_2D_OK.csv'



zscoreCut = 3.09

test = 'Mann-Whitney'
testB = 'Mann-Whitney'

#Apply stat correction ? 'bonferroni' or None
correction = 'bonferroni'

bins = 20

#Do we save the data ?
saveData = False

#Do we save the figure ?
saveFig = False

#Do we perform the stats ?
statsToDo = True

fig, ax = plt.subplots(1,4, figsize=(10,9))
ax[0].set_ylabel('Amplitude (pA+/-SEM)')
ax[0].set_title('Avg synaptic Amplitude')
ax[1].set_ylabel('Total Amplitude (pA+/-SEM)')
ax[1].set_title('Total synaptic Amplitude')
ax[2].set_ylabel('Zscore > {} (%+/-SEM)'.format(zscoreCut))
ax[2].set_title('Proportion of active sites')
ax[3].set_ylabel('Var Coefficient')
ax[3].set_title('Variation coefficient')
plt.tight_layout()

fig2,axx = plt.subplots(1,1)
plt.tight_layout()

fig3,axxx = plt.subplots(1,1)

fig4, box = plt.subplots(1,5,figsize=(12,5))


for i in range(3):
    ax[i].set_xticks(range(len(groups)))
    ax[i].set_xticklabels(groups)
    
def SEM(data):
    return np.nanstd(data)/np.sqrt(len(data))

#FOR THE WHOLE MAP

avg, total, prop, var = [],[],[],[]

MEASURES = []
DISTRIBUTIONS = []

for group,index in zip(groups,range(len(groups))):
    print('Group = {}'.format(group))
    
    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir)]
        
    
    AVG_MEASURE, TOTAL_MEASURE, AVG_PROP, VAR_INDEX = [],[],[],[]
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
        avgMeasure = round(np.nanmean(sitesToKeep),3)
        proportionActiveSite = round(len(sitesToKeep)/len(cleanedMeasures)*100,3)
        totalMeasure = round(np.nansum(sitesToKeep),3)
        
        varIndex = np.nanstd(sitesToKeep) / np.nanmean(sitesToKeep)

        
        #Store distribution
        if idx == 0:
            distributions = sitesToKeep
        else:
            distributions= np.concatenate((distributions,sitesToKeep))
        
        #Append in list
        AVG_MEASURE.append(avgMeasure)
        TOTAL_MEASURE.append(totalMeasure)
        AVG_PROP.append(proportionActiveSite)
        VAR_INDEX.append(varIndex)
        
    
            
            
            
    #Scatter plot             
    ax[0].scatter(np.ones(len(AVG_MEASURE))*index,AVG_MEASURE,color=colors[index])
    ax[0].bar(index,np.nanmean(AVG_MEASURE),alpha=0.5,color=colors[index])
    ax[0].errorbar(index,np.nanmean(AVG_MEASURE),yerr=SEM(AVG_MEASURE),color=colors[index])
    
    avg.append(AVG_MEASURE)
    
    print ('Avg Amplitude per mice at {} = {} +/- {}'.format(group, round(np.nanmean(AVG_MEASURE),3),round(np.nanstd(AVG_MEASURE),3)))

    ax[1].scatter(np.ones(len(TOTAL_MEASURE))*index,TOTAL_MEASURE,color=colors[index])
    ax[1].bar(index,np.nanmean(TOTAL_MEASURE),alpha=0.5,color=colors[index])        
    ax[1].errorbar(index,np.nanmean(TOTAL_MEASURE),yerr=SEM(TOTAL_MEASURE),color=colors[index])
    
    total.append(TOTAL_MEASURE)
    
    print ('Total currents per mice at {} = {} +/- {}'.format(group, round(np.nanmean(TOTAL_MEASURE),3),round(np.nanstd(TOTAL_MEASURE),3)))
    
    ax[2].scatter(np.ones(len(AVG_PROP))*index,AVG_PROP,color=colors[index])
    ax[2].bar(index,np.nanmean(AVG_PROP),alpha=0.5,color=colors[index])
    ax[2].errorbar(index,np.nanmean(AVG_PROP),yerr=SEM(AVG_PROP),color=colors[index])
    
    prop.append(AVG_PROP)
    
    print ('Avg activation per mice at {} = {} +/- {}'.format(group, round(np.nanmean(AVG_PROP),3),round(np.nanstd(AVG_PROP),3)))
    
    ax[3].scatter(np.ones(len(VAR_INDEX))*index,VAR_INDEX,color=colors[index])
    ax[3].bar(index,np.nanmean(VAR_INDEX),alpha=0.5,color=colors[index])
    ax[3].errorbar(index,np.nanmean(VAR_INDEX),yerr=SEM(VAR_INDEX),color=colors[index])
    
    var.append(VAR_INDEX)
    
    axx.hist(distributions,bins=bins,density=True,color=colors[index],label=group,histtype='step')    
    axx.set_ylabel('Norm. Occurence')
    axx.legend(loc='best')
    
    DISTRIBUTIONS.append(distributions)
    
    #Avg values: 
    print ('Avg Amplitude at {} = {} +/- {}'.format(group, round(np.nanmean(distributions),3),round(np.nanstd(distributions),3)))
    
    axxx.hist(distributions,bins=bins,density=True,cumulative=True,color=colors[index],
              histtype='step',label=group)
    
    axxx.set_ylabel('Norm. count')
    axxx.set_xlabel('Amplitude (pA)')
    axxx.legend(loc='best')
    
    #Create a dataframe with all the distributions  
    
    indexList = ['Map#','Condition','Average amplitudes (pA)','Total amplitudes (pA)','Propotion (%)','Var Coeff']

    groupDf = pd.DataFrame(np.vstack((listOfExperiments,
                                      [group for x in range(len(listOfExperiments))],
                                      AVG_MEASURE,
                                      TOTAL_MEASURE,
                                      AVG_PROP,
                                      VAR_INDEX)),index=indexList).transpose()
    
    #Pass numeric columns to float
    groupDf['Average amplitudes (pA)'] = groupDf['Average amplitudes (pA)'].astype(float)
    groupDf['Total amplitudes (pA)'] = groupDf['Total amplitudes (pA)'].astype(float)
    groupDf['Propotion (%)'] = groupDf['Propotion (%)'].astype(float)
    groupDf['Var Coeff'] = groupDf['Var Coeff'].astype(float)
    
    MEASURES.append(groupDf)
    


#Create global dataframe for stats and shit
allData = pd.concat(MEASURES)

distributionsDf = pd.DataFrame(DISTRIBUTIONS, index=groups).T


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

sn.boxplot(x=allData["Condition"],y=allData["Var Coeff"],ax=box[3],palette=pal)
sn.swarmplot(x=allData["Condition"],y=allData["Var Coeff"],ax=box[3],color='0.25')

add_stat_annotation(box[3], x=allData["Condition"], y=allData["Var Coeff"], order=groups,
                    box_pairs=pairs,
                    test=test, text_format='star', loc='inside', verbose=2,comparisons_correction=correction)

sn.boxplot(data=distributionsDf, palette=pal, ax=box[4])
#sn.swarmplot(data=distributionsDf,ax=box[4],color='0.25')

if saveData == True:
    allData.to_excel('{}/Average_Amplitudes.xlsx'.format(saveDir))

if statsToDo == True:

    kruskalAvg = pingouin.kruskal(data = allData, dv='Average amplitudes (pA)',between='Condition')
    print ('Avg Amp')
    print (kruskalAvg)
    box[0].set_title('p(KW)={}'.format(round(kruskalAvg['p-unc'][0],7)))
    
    
    kruskalTotal = pingouin.kruskal(data = allData, dv='Total amplitudes (pA)',between='Condition')
    print ('Total Amp')
    print (kruskalTotal)
    box[1].set_title('p(KW)={}'.format(round(kruskalTotal['p-unc'][0],7)))
    
    kruskalProportion = pingouin.kruskal(data = allData, dv='Propotion (%)',between='Condition')
    print ('Propotion of active sites')
    print (kruskalProportion)
    box[2].set_title('p(KW)={}'.format(round(kruskalProportion['p-unc'][0],7)))
    
    kruskalVar = pingouin.kruskal(data = allData, dv='Var Coeff',between='Condition')
    print ('Variation coefficient')
    print (kruskalVar)
    box[3].set_title('p(KW)={}'.format(round(kruskalVar['p-unc'][0],7)))
    
        
    anovaAvg = pingouin.anova(data = allData, dv='Average amplitudes (pA)',between='Condition')
    print ('Avg Amp')
    print (anovaAvg)
    
    anovaTotal = pingouin.anova(data = allData, dv='Total amplitudes (pA)',between='Condition')
    print ('Total Amp')
    print (anovaTotal)
    
    anovaProportion = pingouin.anova(data = allData, dv='Propotion (%)',between='Condition')
    print ('Propotion of active sites')
    print (anovaProportion)
    
    anovaVar = pingouin.anova(data = allData, dv='Var Coeff',between='Condition')
    print ('Variation coefficient')
    print (anovaVar)
    
    
    
    #STATS
    
    print('----------Stats on amplitude distributions-----------')
    
    allDistKruskal = stats.kruskal(DISTRIBUTIONS[0],
                                   DISTRIBUTIONS[1],
                                   DISTRIBUTIONS[2],
                                   DISTRIBUTIONS[3])
    
    for i in range(len(groups)): 
        
        firstSample = DISTRIBUTIONS[i]
        
        shapiro = stats.shapiro(firstSample)

        print ('Shapiro pvalue on {} = {}'.format(groups[i],shapiro[1]))        
        
        for j in range(len(groups)): 
            
            print('')
            
            secondSample = DISTRIBUTIONS[j]

            ksTest = stats.ks_2samp(firstSample,secondSample,alternative='two-sided',mode='auto')
            
            print ('KS test {} vs {} p value = {}'.format(groups[i],groups[j],ksTest[1]))
            
            t_test = stats.ttest_ind(firstSample,secondSample)

            print ('t test ind test {} vs {} p value = {}'.format(groups[i],groups[j],t_test[1]))
            
            mwu = stats.mannwhitneyu(firstSample,secondSample)
            
            print ('MWU test {} vs {} U value = {}'.format(groups[i],groups[j],mwu[0]))
            print ('MWU test {} vs {} p value = {}'.format(groups[i],groups[j],mwu[1]))
    
    
plt.tight_layout()


if saveFig == True:
    fig.savefig('{}/Amplitudes.pdf'.format(saveDir))
    fig.savefig('{}/Amplitudes.png'.format(saveDir))
    
    fig2.savefig('{}/Amplitudes_Distributions.pdf'.format(saveDir))
    fig2.savefig('{}/Amplitudes_Distributions.png'.format(saveDir))
    
    fig3.savefig('{}/Cumulative_Amplitudes_Distributions.pdf'.format(saveDir))
    fig3.savefig('{}/Cumulative_Amplitudes_Distributions.png'.format(saveDir))
    
    fig4.savefig('{}/BoxPlots.pdf'.format(saveDir))
    fig4.savefig('{}/BoxPlots.png'.format(saveDir))