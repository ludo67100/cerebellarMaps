# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:05:48 2020

@author: ludov
"""

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as stat

file = 'E:/000_PAPER/CatWalk/CATWALK_PROFILES.xlsx'

dataset = pd.read_excel(file,header=0,index_col=0,sheet_name='All')

fig, ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(18,9))
ax[0].set_title('CUFF individual trials') ; ax[0].set_ylabel('Log Balance Index')
ax[1].set_title('SHAM individual trials')


days = dataset.columns[:-1]

for animal in dataset.index[:20]:

    profile = dataset.loc[animal,'baseline':'post_op_33'].values
    
    condition = dataset.loc[animal,'Condition']
    
    if condition == 'CUFF':
        ax[0].plot(profile,label=animal)

    else:
        ax[1].plot(profile,label=animal)

    
for i in range(2):
    ax[i].legend(loc='best')
    ax[i].plot(np.zeros(len(profile)),color='black',linestyle='--')
    ax[i].set_xticks(range(len(days)))
    ax[i].set_xticklabels(days,rotation='45')
    ax[i].set_xlabel('Day#')
      
      
      
# Details : AII and AIII cuff animals seems to describe a 3rd profile
# Do average without these 2 
    
fig2, axx = plt.subplots(1,2,sharex=True,sharey=True,figsize=(18,9))

allCuff = pd.read_excel(file,header=0,index_col=0,sheet_name='AllCuff').values
allCuffAvg = np.mean(allCuff,axis=0)
allCuffSem = stat.sem(allCuff,axis=0)

axx[0].set_title('All individuals')
axx[0].set_ylabel('Avg Log Balance Index +/-SEM')

axx[0].plot(allCuffAvg,color='orange')
axx[0].fill_between(range(len(allCuffSem)),allCuffAvg+allCuffSem,allCuffAvg-allCuffSem,color='orange',alpha=0.2,label='Cuff')

allSham = pd.read_excel(file,header=0,index_col=0,sheet_name='AllSham').values
allShamAvg = np.mean(allSham,axis=0)
allShamSem = stat.sem(allSham,axis=0)


sortedCuff = pd.read_excel(file,header=0,index_col=0,sheet_name='SortedCuff').values
sortedCuffAvg = np.mean(sortedCuff,axis=0)
sortedCuffSem = stat.sem(sortedCuff,axis=0)

axx[1].set_title('AII and AIII excluded')

axx[1].plot(sortedCuffAvg,color='orange')
axx[1].fill_between(range(len(sortedCuffSem)),sortedCuffAvg+sortedCuffSem,sortedCuffAvg-sortedCuffSem,color='orange',alpha=0.2,label='Cuff')


for i in range(2):
    axx[i].plot(allShamAvg,color='black')
    axx[i].fill_between(range(len(allShamSem)),allShamAvg+allShamSem,allShamAvg-allShamSem,color='black',alpha=0.2,label='Sham')
    
    axx[i].legend(loc='best')
    axx[i].plot(np.zeros(len(profile)),color='black',linestyle='--')
    axx[i].set_xticks(range(len(days)))
    axx[i].set_xticklabels(days,rotation='45')
    axx[i].set_xlabel('Day#')
