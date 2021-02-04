# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:38:33 2020

@author: ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
import pandas as pd
import seaborn as sn
sn.set_style("whitegrid")

from scipy import stats
from scipy.stats import median_absolute_deviation as MAD
import math
import pingouin as pg 
import sys 


#Define the groups
groups = ['P9P10','P12P13','P14P18','P30P40']
colors = ['lightskyblue','skyblue','deepskyblue','royalblue']
pal = ['lightskyblue','skyblue','deepskyblue','royalblue']


#Input folder
dataSource = 'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS'


#deepValues = ['Upper','Up-Mid','Up-Low','Lower']
deepValues = ['Upper','Lower']

#Savedir 
saveDir =  'E:/000_PAPER/Development/Amplitude_Analysis/11_Deepness_Analysis'

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'
positionFileType = 'Positions_cp_centered_OK.csv' 

zscoreCut = 3.09

#Do we save the data ?
saveData = True

#Do we save the figure ?
saveFig = True

#Round factor (decimals)
RF = 2

#Redirect console output to a txt file instead of console 
sys.stdout = open('{}/Statistics.txt'.format(saveDir),'w')

    
def SEM(data):
    return np.nanstd(data)/np.sqrt(len(data))

#def activation_by_line(AmpMap, Zscore, threshold=3.09):
#    
#    temp = []
#    
#    
#    for line in range(AmpMap.shape[0]):
#        
#        significant_sites = [x for x in Zscore[line,:].ravel() if math.isnan(x)==False and x >= threshold]
#        total_sites = [x for x in AmpMap[line,:].ravel() if math.isnan(x)==False]
#        
#        perc = len(significant_sites) / len(total_sites) * 100.
#        
#        temp.append(perc)
#        
#    
#    return temp


def activation_by_line(AmpMap, Zscore, threshold=3.09):
    
    temp = []
    
    lines = [[0,2],[2,4]]
    
    for line in lines:
        
        significant_sites = [x for x in Zscore[line[0]:line[1],:].ravel() if math.isnan(x)==False and x >= threshold]
        total_sites = [x for x in AmpMap[line[0]:line[1],:].ravel() if math.isnan(x)==False]
        
        perc = len(significant_sites) / len(total_sites) * 100.
        
        temp.append(perc)
        
    
    return temp



#FOR THE WHOLE MAP

globalDfAvg = []
globalDfStd = []
globalDf = []

boxplot, ax = plt.subplots(1,len(groups), sharex=True, sharey=True, figsize=(16,3))
ax[0].set_ylabel('GC section #')


for group,index in zip(groups,range(len(groups))):
    print ('')
    print('Group = {}'.format(group))
    
    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir)]
        
    listOfLineActivation = []    
    computedExperiment = 0

    for manip,idx in zip(listOfExperiments,range(len(listOfExperiments))):
                
        #Get amplitudes or charges
        measures = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,fileType),delimiter=','))
        #Get corresponding zscores
        zscores = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,zscoreFileType),delimiter=','))
        #Get positions 
        positions = np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,positionFileType),delimiter=',')
        
        listOfLineActivation.append(activation_by_line(measures,zscores))
    
    activationDf = pd.DataFrame(listOfLineActivation,index=listOfExperiments,
                                columns=deepValues)
    
    sn.boxplot(data=activationDf, orient='h',ax=ax[index],palette=[pal[index] for x in range(len(deepValues))])
    ax[index].set_title(group)
    ax[index].set_xlabel('% active sites')

    
    globalDfAvg.append(activationDf.mean())
    globalDfStd.append(activationDf.std())
    
    #Create a copy of activationDf with Condition column
    activationDfcopy = activationDf.copy(deep=True)
    activationDfcopy['condition'] = [groups[index] for x in range(activationDfcopy.shape[0])]
    globalDf.append(activationDfcopy)
    
    
    #Intra group statistics
    print ('---- Intra Group Statistics ----')
    
    print ('Avg activation (%) in {}: {}+/-{}'.format(deepValues[0],
           round(activationDf[deepValues[0]].mean(),RF), 
           round(activationDf[deepValues[0]].std(),RF)))
    
    print ('Avg activation (%) in {}: {}+/-{}'.format(deepValues[1],
           round(activationDf[deepValues[1]].mean(),RF), 
           round(activationDf[deepValues[1]].std(),RF)))
    
    print ('Med activation (%) in {}: {}+/-{}'.format(deepValues[0],
           round(activationDf[deepValues[0]].median(),RF), 
           round(MAD(activationDf[deepValues[0]].values),RF)))
    
    print ('Med activation (%) in {}: {}+/-{}'.format(deepValues[1],
           round(activationDf[deepValues[1]].median(),RF), 
           round(MAD(activationDf[deepValues[1]].values),RF)))
    
    kwTest = stats.kruskal(activationDf[deepValues[0]],
                           activationDf[deepValues[1]])
    
    print ('Kruskal Wallis test H value = {}'.format(kwTest[0]))
    print ('Kruskal Wallis test p value = {}'.format(kwTest[1]))
    
    if kwTest[1] <=0.05:
        
        print('Post Hoc MWU test')
    
        #Do pairwise MWU tests: 
        
        for indexA,colA in zip(range(activationDf.shape[1]),activationDf.columns):
            
            refCol = activationDf.iloc[:,indexA]
            
            for indexB, colB in zip(range(activationDf.shape[1]),activationDf.columns):
                
                compareCol = activationDf.iloc[:,indexB]
                
                mwu_test = stats.mannwhitneyu(refCol,compareCol,alternative='two-sided')
                if mwu_test[1] <= 0.05/float(activationDf.shape[1]):
                    print ('    {} vs {} U-val={}'.format(colA, colB,mwu_test[0]))
                    print ('    {} vs {} p-val={}'.format(colA, colB,mwu_test[1]))
    
    else:
        print ('KW test failed')
    
    
    
plt.tight_layout()

#Inter group stats 
print ('')
print ('------Inter Group Statistics-------')
globalDf = pd.concat(globalDf,axis=0)

for col in deepValues: 
    
    interKruskal = pg.kruskal(data=globalDf, dv=col,between='condition')
    print ('Range:{}'.format(col))
    print (interKruskal)
    print('')
    
    #post hoc MWU
    for groupA in groups:
        serieA = globalDf[col].loc[(globalDf['condition']==groupA)]
        
        for groupB in groups: 
            serieB = globalDf[col].loc[(globalDf['condition']==groupB)]
            
            interPostMwu = stats.mannwhitneyu(serieA, serieB,alternative='two-sided')
            
            if interPostMwu[1] >= 0.05:
                mwuResult = 'FAIL'
            else:
                mwuResult = 'PASS'
            
            print ('post hoc MWU {} vs {}; U={}     '.format(groupA, groupB, interPostMwu[0]))
            print ('post hoc MWU {} vs {}; p={}     {}'.format(groupA, groupB, interPostMwu[1],mwuResult))
            

    



#KSavings
if saveFig == True: 
    boxplot.savefig('{}/GC_Layer_Activation.pdf'.format(saveDir))
    boxplot.savefig('{}/GC_Layer_Activation.png'.format(saveDir))    
    
if saveData == True:
    AvgDf = pd.concat(globalDfAvg,axis=1)
    AvgDf.columns=groups
    
    
    StdDf = pd.concat(globalDfStd,axis=1)
    StdDf.columns = groups
    
    writer = pd.ExcelWriter('{}/GC_Deepness_Actication.xlsx'.format(saveDir))
    
    AvgDf.to_excel(writer, sheet_name='AVG ACTIVATION')
    StdDf.to_excel(writer, sheet_name='STD')
    
    writer.close()
    
    
    globalDf.to_excel('{}/Deepness_values.xlsx'.format(saveDir))



    
    
        
        
        