# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 18:25:16 2020

@author: ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd 
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt 
import numpy as np 
from scipy.stats import median_absolute_deviation as MAD


#t-SNE params to play with 
perp = 10
metric = 'euclidean'


#I/O path
file = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/DataSource_Spaeth_Bahuguna_et_al.xlsx'
savedir = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/REVISION CODE/00_FINAL_CODE/Figure5/5A'
colors =['skyblue','orange','purple','lightcoral','0.5','green','limegreen']


#Get data
dataset = pd.read_excel(file,index_col=0,header=1,sheet_name='Fig5a-input')

data = dataset.drop(columns='Condition')

#Freeze the seed to reproduce result
seed = 656610
np.random.seed(seed)

#Init tsne
tsne = TSNE(n_components=2, verbose=1, perplexity=perp, n_iter=3000,init='random',metric=metric,learning_rate=150.0,early_exaggeration=24.)

#Fit with data
fitTSNE = pd.DataFrame(tsne.fit_transform(data.values),columns=['X t-SNE', 'Y t-SNE'])

fitTSNE['Condition'] = dataset['Condition'].values

#Init plot
plt.figure(figsize=(6,6)) ; plt.xlabel('x t-SNE'); plt.ylabel('y t-SNE')
plt.title('t-SNE on synaptic params (perplexity={} / metric={})'.format(perp,metric))

#Plot for each condition
for cond, i in zip(np.unique(dataset['Condition'].values),range(len(np.unique(dataset['Condition'].values)))):
    
    
    #Slice data according to condition
    subDf = fitTSNE.loc[fitTSNE['Condition']==cond]
    plt.scatter(subDf['X t-SNE'],subDf['Y t-SNE'],alpha=0.2,color=colors[i])
    
    #Compute cloud kernel (median)
    kernel = [np.median(subDf['X t-SNE'].values), np.median(subDf['Y t-SNE'].values)]
    quartilesX = list(np.abs(kernel[0]-np.array([np.percentile(subDf['X t-SNE'].values,25),np.percentile(subDf['X t-SNE'].values,75)])))
    quartilesY = list(np.abs(kernel[1]-np.array([np.percentile(subDf['Y t-SNE'].values,25),np.percentile(subDf['Y t-SNE'].values,75)])))
    plt.scatter(kernel[0],kernel[1],color=colors[i],label=cond,s=70)

    #Compute errorbars (MAD)
    MADs = MAD(subDf.drop(columns='Condition').values)
    plt.errorbar(kernel[0],kernel[1],xerr=np.array([[abs(quartilesX[0]),abs(quartilesX[1])]]).T,
                 yerr=np.array([[abs(quartilesY[0]),abs(quartilesY[1])]]).T,color=colors[i])
    
plt.legend(loc='best')
plt.savefig('{}/t-SNE_synaptic_params_perp{}_metric-{}.pdf'.format(savedir, perp, metric))