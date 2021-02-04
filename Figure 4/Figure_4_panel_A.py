# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 19:29:22 2020

@author: Ludovic.spaeth
"""


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

import numpy as np
import matplotlib.pyplot as plt
import os 
import pandas as pd
import seaborn as sn
from statannot import add_stat_annotation
from scipy import stats
import sys
from scipy.stats import median_absolute_deviation as MAD




groups = ['WT','ENR1','ENR2','LC','LS','EC','ES']
colors = ['skyblue','limegreen','green','lightcoral','black','orange','purple']
pal = ['skyblue','limegreen','green','lightcoral','grey','orange','purple']

#Pairs to analyse for stats
pairs= [("WT", "ENR1"), ("WT", "ENR2"),("WT", "LC"), ("WT", "LS"), ("WT", "EC"),("WT", "ES")]


#Input folder
dataSource = 'E:/000_PAPER/Amplitude_Analysis/Sigma'

#Savedir  
saveDir =  'E:/000_PAPER/Amplitude_Analysis/Sigma/09_AMPLITUDE_DISTRIBUTIONS'


#Redirect console output to a txt file instead of console 
sys.stdout = open('{}/Syn_Amp_Statistics.txt'.format(saveDir),'w')

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'


zscoreCut = 3.09

test = 'Mann-Whitney'
testB = 'Mann-Whitney'

#Apply stat correction ? 'bonferroni' or None
correction = None

bins = 200

#Do we save the data ?
saveData = False

#Do we save the figure ?
saveFig = False

#Do we perform the stats ?
statsToDo = True

#Round factor (decimals)
RF =2

figure, ax = plt.subplots(1,1)

#Subplot for distributions 
curveplot, curves = plt.subplots(1,len(groups),sharex=True, sharey=True, figsize=(10,2))
    
def SEM(data):
    return np.nanstd(data)/np.sqrt(len(data))

distributions = []

#FOR THE WHOLE MAP
for group,index in zip(groups,range(len(groups))):
    print('')
    print('Group = {}'.format(group))
    
    currents = []
    
    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir)]
        
    #Get currents for each manip
    for manip,idx in zip(listOfExperiments,range(len(listOfExperiments))):
                
        #Get amplitudes or charges
        measures = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,fileType),delimiter=',').ravel())
        #Get corresponding zscores
        zscores = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,zscoreFileType),delimiter=',').ravel())

        #Extract significant measures on the whole map 
        sitesToKeep = [x for (x,y) in zip(measures,zscores) if y >= zscoreCut]
        
        for i in sitesToKeep:
            currents.append(i)
    
    #Average values
    distributions.append(currents)
    
    avg = np.nanmean(currents)
    std = np.nanstd(currents)
    
    med = np.nanmedian(currents)
    mad = MAD(currents)
    
    print ('MEAN {} -> {} +/- {} pA'.format(group, round(avg,RF), round(std,RF)))
    print ('MEDIAN {} -> {} +/- {} pA'.format(group, round(med,RF), round(mad,RF)))
    print ('Min. current: {} pA'.format(min(currents)))
    print ('Max. current: {} pA'.format(max(currents)))
    
    sn.distplot(currents, kde=True, color=colors[index],ax=ax)
    sn.distplot(currents, hist=True, kde=True, color=colors[index], ax=curves[index])
    
#Do a dataframe
data = pd.DataFrame(distributions, index=groups).T

boxplot, box = plt.subplots(1,1)
box.set_ylabel('Amplitude (pA)')

sn.boxplot(data=data, palette=pal, ax=box)


#Do a small plot with inset on average values for each values 
groups = ['WT','ENR1','ENR2','LC','LS','EC','ES']

little_plot, axplot = plt.subplots(1,1, figsize=(5,3))
axplot.set_ylabel('Avg Amplitude (pA)')

x_axis_indexes = [0,1,2,2,2,1,1]

for i in range(len(x_axis_indexes)): 
    
    axplot.scatter(x_axis_indexes[i],np.nanmean(distributions[i]),color=colors[i],label=groups[i])
    
    axplot.errorbar(x_axis_indexes[i],np.nanmean(distributions[i]),color=colors[i],yerr=SEM(distributions[i])) 

axplot.legend(loc='best',ncol=3)
    
#----------------------------------------------
ctrlEnrKruskal = stats.kruskal(distributions[0],
                               distributions[1],
                               distributions[2])
print('')
print ('Non parametric stats on CTRL/ENR - CTRL/CUFF - CTRL/SHAM')
print ('CTRL-ENR1-ENR2 Kruskal Wallis test H value = {}'.format(ctrlEnrKruskal[0]))
print ('CTRL-ENR1-ENR2 Kruskal Wallis test p value = {}'.format(ctrlEnrKruskal[1]))

if ctrlEnrKruskal[1] <= 0.05: 
    mwuCtrlEnr1 = stats.mannwhitneyu(distributions[0],distributions[1],alternative='two-sided')
    mwuCtrlEnr2 = stats.mannwhitneyu(distributions[1],distributions[2],alternative='two-sided')

    print ('CTRL vs ENR1 mwu test U value = {}'.format(mwuCtrlEnr1[0]))
    print ('CTRL vs ENR1 mwu test p value = {}'.format(mwuCtrlEnr1[1]))
    print ('ENR1 vs ENR2 mwu test U value = {}'.format(mwuCtrlEnr2[0]))
    print ('ENR1 vs ENR2 mwu test p value = {}'.format(mwuCtrlEnr2[1]))

else:
    print('KW test failed')

#-----------------------------------------------
ctrlShamKruskal = stats.kruskal(distributions[0],
                                distributions[4],
                                distributions[6])
print('')
print ('CTRL-ES-LS Kruskal Wallis test H value = {}'.format(ctrlShamKruskal[0]))
print ('CTRL-ES-LS Kruskal Wallis test p value = {}'.format(ctrlShamKruskal[1]))

if ctrlShamKruskal[1] <= 0.05: 
    mwuCtrlEs = stats.mannwhitneyu(distributions[0],distributions[6],alternative='two-sided')
    mwuEsLS = stats.mannwhitneyu(distributions[4],distributions[6],alternative='two-sided')
    
    print ('CTRL vs ES mwu test U value = {}'.format(mwuCtrlEs[0]))
    print ('CTRL vs ES mwu test p value = {}'.format(mwuCtrlEs[1]))
    print ('ES vs LS mwu test U value = {}'.format(mwuEsLS[0]))
    print ('ES vs LS mwu test p value = {}'.format(mwuEsLS[1]))
    
else:
    print('KW test failed')
    
#-----------------------------------------------
ctrlCuffKruskal = stats.kruskal(distributions[0],
                                distributions[3],
                                distributions[5])
print('')
print ('CTRL-EC-LC Kruskal Wallis test H value = {}'.format(ctrlCuffKruskal[0]))
print ('CTRL-EC-LC Kruskal Wallis test p value = {}'.format(ctrlCuffKruskal[1]))

if ctrlCuffKruskal[1] <= 0.05: 
    mwuCtrlEc = stats.mannwhitneyu(distributions[0],distributions[5],alternative='two-sided')
    mwuEcLc = stats.mannwhitneyu(distributions[3],distributions[5],alternative='two-sided')

    print ('CTRL vs EC mwu test U value = {}'.format(mwuCtrlEc[0]))  
    print ('CTRL vs EC mwu test p value = {}'.format(mwuCtrlEc[1]))
    print ('EC vs LC mwu test U value = {}'.format(mwuEcLc[0]))
    print ('EC vs LC mwu test p value = {}'.format(mwuEcLc[1]))
    
else:
    print('KW test failed')


#Statistics 
#-----------------Non parametric-----------------------------------------------
kruskal = stats.kruskal(distributions[0],
                        distributions[1],
                        distributions[2],
                        distributions[3],
                        distributions[4],
                        distributions[5])

print('')
print ('Non parametric stats on every conditions')
print ('Kruskal-Wallis H value={}'.format(kruskal[0]))
print ('Kruskal-Wallis p value={}'.format(kruskal[1]))

#Post hoc mwu tests
for i in range(len(groups)):
    ref = distributions[0] #WT group
    
    mwu = stats.mannwhitneyu(ref, distributions[i], alternative='two-sided')
    
    print ('MWU {} vs {} Uvalue={}'.format(groups[0], groups[i], mwu[0]))    
    print ('MWU {} vs {} pvalue={}'.format(groups[0], groups[i], mwu[1]))

#---------------Parametric-----------------------------------------------------
anova = stats.f_oneway(distributions[0],
                       distributions[1],
                       distributions[2],
                       distributions[3],
                       distributions[4],
                       distributions[5])

print ('Parametric stats')

print ('OneWay Anova p value={}'.format(anova[1]))

#Post hoc t tests
for i in range(len(groups)):
    ref = distributions[0] #WT group
    
    ttest = stats.ttest_ind(ref, distributions[i], nan_policy='omit')
    print ('ind. T-test {} vs {} pvalue={}'.format(groups[0], groups[i], ttest[1]))
    
    welch = stats.ttest_ind(ref, distributions[i], equal_var=False, nan_policy='omit')
    print ('Welch T-test {} vs {} pvalue={}'.format(groups[0], groups[i], welch[1]))
    
    
    
    
