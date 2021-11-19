# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 16:12:24 2020

@author: ludov
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sn 
from scipy import stats

from scipy.optimize import curve_fit


file = 'C:/Users/klab/Downloads/DataSource_Spaeth_Bahuguna_et_al (1).xlsx'


if __name__ == '__main__':

    dataset = pd.read_excel(file,header=1,index_col=0,na_values='n/a',sheet_name='Figure 3a').iloc[:13]
    
    fig, ax = plt.subplots(1,1)
    plt.title('ENR individual plots')
    
    ax.set_ylabel('Distance per hour (m)')
    
    profiles = []
    
    
    for animal in dataset.index: 
        
        profile = dataset.loc[animal,'Day1':'Day19']
        
        plt.plot(profile,color='0.5',alpha=0.3)
        
        profiles.append(profile.dropna().values)
        
    
    Avg = dataset.mean(axis=0).dropna()
    Std = dataset.std(axis=0).dropna()
    
    ax.plot(Avg, linewidth=2,color='green')
    ax.fill_between(np.arange(0,len(Avg),1),Avg+Std, Avg-Std, alpha=0.2,color='green')
    
    #Test day 1 vs 19 with paired t-test
    ttest = stats.ttest_rel(dataset.loc[dataset['Condition']=='long training']['Day1'].values,
                            dataset.loc[dataset['Condition']=='long training']['Day19'].values)

    print('paired T-test day 1 - 19')
    print(ttest)    
    

    X,Y = np.asarray([i for i in range(1,20)]),Avg.values
    
    def func_mono_exp(x, a, b, c, d):
        return a * np.exp(-(x-b)/c) + d
    
    #do mono exp fit
    popt,pcov = curve_fit(func_mono_exp,X,Y)

    #Plot fit
    newX = np.unique(X)
    #ax.plot(newX,expFit[0]*np.log(newX)+expFit[1], label='fit [a*log(x)+b]', linewidth=2)
    
    projX = np.arange(0,20,1)
    ax.plot(func_mono_exp(projX,popt[0],popt[1],popt[2],popt[3]),label='fit [a*e^(-(x-b)/c)]',linewidth=2)
    
    #Determine horizontal asymptote of exp function 
    y_asymptote = np.nanmean(func_mono_exp(np.arange(0,1000000,1),popt[0],popt[1],popt[2],popt[3]))
    
    #plot asymptote 
    ax.plot(projX,np.ones(len(projX))*y_asymptote,color='black',linestyle='--',label='lim (0->+inf)={}'.format(round(y_asymptote,2)))
    
    ax.legend(loc='best')
    
    #Check goodness of fit with exp 
    ksStat = stats.ks_2samp(func_mono_exp(projX,popt[0],popt[1],popt[2],popt[3]),Y)
    
    print ('KS test between exponential fit and data, p={}'.format(ksStat[1]))
    
    if ksStat[1] > 0.05:
        print('Exp. fit and Data cannot be discriminated : perfect fit')
    else:
        print('Exp. fit and Data are different distributions : weak fit')


#    
#    #Plot extended fit for asymptot
#    assX = np.arange(0,50,1)
#    ax.plot(assX,expFit[0]*np.log(assX)+expFit[1], label='extended fit', linewidth=2,linestyle='--')
#    
#    ax.legend(loc='best')
#    
#    
#    a = Fit_single_trace(Y,X,0,10)
        
    #plt.legend(loc='best')
    
    
#    plt.plot(newX,expFit[0]*np.log(newX)+expFit[1], label='fit (a*log(x)+b)', linewidth=2)
#    
#    newY = expFit[0]*np.log(newX)+expFit[1]
#    
#    plt.plot(newX, np.exp((newY-expFit[1])/expFit[0]), label='fit (e^y-b/a)')
#    
#    plt.legend(loc='best')


#Plot the behavioral features 

features = pd.read_excel(file,header=17,index_col=0,na_values='n/a',sheet_name='Figure 3a')

fig2, axx  = plt.subplots(1,3)


featureLabels = ['Intercept','Slope','Total distance']

for i in range(len(featureLabels)):
    
    sn.boxplot(x='Condition',y=featureLabels[i],data=features,ax=axx[i])
    sn.stripplot(x='Condition',y=featureLabels[i],data=features,ax=axx[i],color='black')
    
    x = features.loc[features['Condition']=='short training'][featureLabels[i]].values
    y = features.loc[features['Condition']=='long training'][featureLabels[i]].values
    
    mwu = stats.mannwhitneyu(x,y, alternative = 'two-sided')
    print('')
    print ('{} - short training avg +/- sd (n={}) = {:.2f} +/- {:.2f}'.format(featureLabels[i], len(x), np.nanmean(x), np.nanstd(x)))
    print ('{} - long training avg +/- sd (n={}) = {:.2f} +/- {:.2f}'.format(featureLabels[i], len(y), np.nanmean(y), np.nanstd(y)))
    print ('{} : short vs long training'.format(featureLabels[i]))
    print ('Two-tailed MWU, U={} ; p={}'.format(mwu[0],mwu[1]))
    
plt.tight_layout()



    

