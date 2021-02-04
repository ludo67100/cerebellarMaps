# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 18:54:00 2020

@author: ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sn 
import pingouin as pg 
from scipy import stats
import sys


groups = ['WT','ENR1','ENR2','LC','LS','EC','ES']
colors = ['skyblue','limegreen','green','lightcoral','0.3','orange','purple']




file = 'E:/000_PAPER/FIGURES/Synaptic/LocvsDist/Activation_Loc_vs_Dist.xlsx'

saveDir = 'E:/000_PAPER/FIGURES/Synaptic/LocvsDist/Loc vs Loc and Dist vs Dist Activation'

sys.stdout = open('{}/Active_sites_LocVSLoc_DistVSDist_Statistics.txt'.format(saveDir),'w')

fig, ax = plt.subplots(1,2,sharex=True,sharey=True)
plt.suptitle('Activation - Local vs Distal')

df = pd.read_excel(file, header=0, index_col=0)

df.index=df['condition']

df = df.loc[groups]


for cat,indx in zip(['Local','Distal'],range(2)): 
    
    print ('')
    print ('{} side'.format(cat))
    
    sn.boxplot(x='condition',y=cat,data=df,ax=ax[indx],palette=colors)
    sn.swarmplot(x='condition',y=cat,data=df,ax=ax[indx],color='0.5')
    
    print(pg.kruskal(data=df, dv=cat,between='condition'))
    
    print('')
    
    
#    kw = stats.kruskal(df.loc[groups[0],cat].values,
#                       df.loc[groups[1],cat].values,
#                       df.loc[groups[2],cat].values,
#                       df.loc[groups[3],cat].values,
#                       df.loc[groups[4],cat].values,
#                       df.loc[groups[5],cat].values,
#                       df.loc[groups[6],cat].values)
    
#    print (kw)
    
    for group in groups: 
        
        refGroup = df.loc['WT',cat]
        
        compGroup = df.loc[group,cat]
        
        mwu = stats.mannwhitneyu(refGroup,compGroup)
        
        print ('MWU WT vs {} U-val={}'.format(group,mwu[0]))
        print ('MWU WT vs {} p-val={}'.format(group,mwu[1]))
        
plt.savefig('{}/Ipsi_Contra.pdf'.format(saveDir))
plt.savefig('{}/Ipsi_Contra.png'.format(saveDir))
        
        

