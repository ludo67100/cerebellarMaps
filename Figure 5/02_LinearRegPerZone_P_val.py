# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:02:43 2020

@author: Ludovic.spaeth
"""


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sn

import pandas as pd 
import numpy as np

from matplotlib import pyplot as plt

from scipy import stats as stat


url = 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max'
savedir= 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max/FIGURES'

conditions = ['WT','ENR1','ENR2','LC','LS','EC','ES']
colors = ['skyblue','limegreen','green','lightcoral','black','orange','purple']
cmaps = ['Blues','YlGn','Greens','Reds','Greys','Oranges','Purples']
rcmaps = ['Blues_r','YlGn_r','Greens_r','Reds_r','Greys_r','Oranges_r','Purples_r']


winsorize = False
cutoff = 0.05 #For winsorize

N = 8
r_min = 0.49999 #Minimum r value for detection as correlation 
zscore = 3.09

savefig = True
savedata = True

vmin = r_min

individual_plot = False 

savefigindividualplots = False

#----------------------------------FUNCTIONs-----------------------------------
#------------------------------------------------------------------------------
def LinReg(x,y,conf=0.95,plot=True):
    '''
    Computes linear regression on 2 arrays 
    
    x,y (1D arrays) = the x & y data 
    
    conf (float) = confidence treshold (alpha = 1-conf)
    
    plot (bool) : if True, will plot the result of the linear reg 
    
    Returns : 
        
        px (array) : x-axis for plot, based on min/max values of x
        
        nom (array) : y-values of linear model (nomial values of y)
        
        lpb, upb (arrays) : lower and upper prediction bands 
        
        r2 = squared correlation coefficient 
    
    '''
    import numpy as np 
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    from scipy import stats
    import uncertainties as unc
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42
    
    #Pip install uncertainties if needed 
    try :
        import uncertainties.unumpy as unp
        
    except : 
        import pip
        pip.main(['install','uncertainties'])
        import uncertainties.unumpy as unp 
    
    
    n = len(y)
    
    def f(x, a, b):
        return a * x + b
    
    popt, pcov = curve_fit(f, x, y)
    
    # retrieve parameter values
    a = popt[0]
    b = popt[1]

    # compute r^2
    r2 = 1.0-(sum((y-f(x,a,b))**2)/((n-1.0)*np.var(y,ddof=1)))
    
    # calculate parameter confidence interval
    a,b = unc.correlated_values(popt, pcov)

    
    # plot data
    if plot == True:
        plt.figure()
        plt.scatter(x, y, s=20, label='Data')
    
    # calculate regression confidence interval
    px = np.linspace(np.min(x),np.max(x),n)
    py = a*px+b
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)
    
    def predband(x, xd, yd, p, func, conf=conf):
        '''
        x = requested points
        xd = x data
        yd = y data
        p = parameters
        func = function name
        '''
        alpha = 1.0 - conf    # significance
        N = xd.size          # data sample size
        var_n = len(p)  # number of parameters
        # Quantile of Student's t distribution for p=(1-alpha/2)
        q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
        # Stdev of an individual measurement
        se = np.sqrt(1. / (N - var_n) * \
                     np.sum((yd - func(xd, *p)) ** 2))
        # Auxiliary definitions
        sx = (x - xd.mean()) ** 2
        sxd = np.sum((xd - xd.mean()) ** 2)
        # Predicted values (best-fit model)
        yp = func(x, *p)
        # Prediction band
        dy = q * se * np.sqrt(1.0+ (1.0/N) + (sx/sxd))
        # Upper & lower prediction bands.
        lpb, upb = yp - dy, yp + dy
        return lpb, upb
    
    lpb, upb = predband(px, x, y, popt, f, conf=conf)
    
    if plot == True:
        # plot the regression
        plt.plot(px, nom, c='orange', label='y=a x + b',linewidth=2)
        
        # uncertainty lines (95% confidence)
        plt.plot(px, nom - 1.96 * std, c='0.5',linestyle='--',\
                 label='95% Confidence Region')
        plt.plot(px, nom + 1.96 * std, c='0.5',linestyle='--')
        # prediction band (95% confidence)
        plt.plot(px, lpb, color='0.5',label='95% Prediction Band',linestyle=':')
        plt.plot(px, upb, color='0.5',linestyle=':')
        plt.ylabel('Y')
        plt.xlabel('X')
        plt.legend(loc='best')
        plt.title('Linear Reg : R$^2$={}'.format(round(r2,2)))
        plt.show()
        
    return px,nom,lpb,upb,r2,std


#-------------------------------------THE SERIOUS SHIT-----------------------------------
#------------------------------------INCLUDING OUTLIERS-------------------------------------------

for cond,color,cmap,rcmap in zip(conditions,colors,cmaps,rcmaps) :
    
    print ('------------------{}---------------------'.format(cond))
    
    try : 
    
        Rs = np.zeros((N,N))
        Ps = np.zeros((N,N))
        
        fig, ax = plt.subplots(N,N, figsize=(16,9))
        fig.suptitle('{} dataset : average amplitude in cluster correlation'.format(cond))
        fig.subplots_adjust(wspace = 0.5,hspace=0.4)
        
        fig2,axx = plt.subplots(1,1,figsize=(N,N))
        
        fig3,axxx = plt.subplots(1,1,figsize=(N,N))
        
        path = '{}/{}_Map_2D_Average_Amp_per_zones_wNANS_{}.xlsx'.format(url,cond,zscore) 
    
        dataset = pd.read_excel(path)
    
        bands = dataset.columns[1:]
        animals = dataset.index
        
        for band_A,y_idx in zip(bands,range(len(bands))):
    
            #Get values from the bands 
            raw_y = np.ravel(dataset[[band_A]].values)
            
            for band_B, x_idx in zip(bands,range(len(bands))):
                print (band_A, ' vs ', band_B)
                
                raw_x = np.ravel(dataset[[band_B]].values)
                
                x,y = [],[]
                
                for i in range(len(raw_y)):
                    if np.isnan(raw_x[i])==False and np.isnan(raw_y[i])==False:
                        x.append(raw_x[i])
                        y.append(raw_y[i])
                        
                if winsorize == True :
                    x = stat.winsorize(np.asarray(x),limits=[cutoff, cutoff])
                    y = stat.winsorize(np.asarray(y),limits=[cutoff, cutoff])
                    
                else: 
                    x = np.asarray(x)
                    y = np.asarray(y)
                    
                
                px,nom,lpb,upb,r2,std = LinReg(x,y,conf=0.95,plot=False)
                
                linearRegStat = stat.linregress(x,y)
                #print (linearRegStat[3])
                

                if r2>=0.99:
                    title_color='0.5'  
                elif np.sqrt(r2) >= r_min:
                    title_color='r'
                else:
                    title_color='black'
                

                
                ax[x_idx,y_idx].set_title('R={}, R$^2$={}'.format(round(np.sqrt(r2),2),round(r2,2)),color=title_color)
                ax[x_idx,y_idx].scatter(x,y,color=color,s=7)
    
                ax[x_idx,y_idx].plot(px,nom,color=color)
                ax[x_idx,y_idx].set_xticks([]);ax[x_idx,y_idx].set_yticks([])
                if x_idx == len(bands)-1:
                    ax[x_idx,y_idx].set_xlabel(band_A)
                if y_idx ==0: 
                    ax[x_idx,y_idx].set_ylabel(band_B)
                
                ax[x_idx,y_idx].plot(px, nom - 1.96 * std, c='0.5',linestyle='--')
                ax[x_idx,y_idx].plot(px, nom + 1.96 * std, c='0.5',linestyle='--')
                ax[x_idx,y_idx].plot(px, lpb, color='0.5',linestyle=':')
                ax[x_idx,y_idx].plot(px, upb, color='0.5',linestyle=':')  
                
                if individual_plot == True : 
                    fig4 = plt.figure()
                    plt.suptitle('R={}, R$^2$={}'.format(round(np.sqrt(r2),2),round(r2,2)),color=title_color)
                    plt.scatter(x,y,color=color,s=7)
        
                    plt.plot(px,nom,color=color)
                    #plt.xticks([]);plt.yticks([])

                    plt.xlabel(band_A)

                    plt.ylabel(band_B)
                    
                    plt.plot(px, nom - 1.96 * std, c='0.5',linestyle='--')
                    plt.plot(px, nom + 1.96 * std, c='0.5',linestyle='--')
                    plt.plot(px, lpb, color='0.5',linestyle=':')
                    plt.plot(px, upb, color='0.5',linestyle=':')   
                    
                    if savefigindividualplots == True:
                        plt.savefig('{}/{}_{}_vs_{}_scatter.pdf'.format(savedir,cond,band_A, band_B))
                        
                    plt.close()
                    

                Ps[x_idx,y_idx] = linearRegStat[3]

                
                if np.sqrt(r2) <= .99:
                    Rs[x_idx,y_idx] = np.sqrt(r2)
                else :
                    Rs[x_idx,y_idx] = 0.0
                    
                
                
                
        axx.imshow(Rs>=r_min,cmap=cmap, vmax=np.max(Rs))
        axx.set_title('{} corr coefficients'.format(cond))
        
        axxx.imshow(Ps,cmap=rcmap, vmax=0.05)
        axxx.set_title('{} P_values'.format(cond))
        
        R_coeffs = []
        R_labels = []
        P_values = []
        
        for (x, y),i in np.ndenumerate(Rs):

            
            if i >= r_min:
                value_color = 'white'
            else:
                value_color = '0.2'
            
            text = axx.text(x,y,round(i,2),ha='center',va='center',color=value_color)
            

            R_coeffs.append(i)
            R_labels.append('{}vs{}'.format(bands[x], bands[y]))
            
        for (x, y),i in np.ndenumerate(Ps):
            text = axxx.text(x,y,round(i,4),ha='center',va='center')
            
            P_values.append(i)
            
        
        
        
        axx.set_xticks(np.arange(0,len(bands),1))
        axx.set_yticks(np.arange(0,len(bands),1))
        axx.set_xticklabels(bands,rotation='vertical')
        axx.set_yticklabels(bands)   
        
        axxx.set_xticks(np.arange(0,len(bands),1))
        axxx.set_yticks(np.arange(0,len(bands),1))
        axxx.set_xticklabels(bands,rotation='vertical')
        axxx.set_yticklabels(bands)   
        
        if savefig == True:
            fig.savefig('{}/{}_linear_regs_zscore{}.pdf'.format(savedir,cond,zscore))
            fig2.savefig('{}/{}_zones_correlations_zscore{}.pdf'.format(savedir,cond,zscore))
            fig3.savefig('{}/{}_pvalues_zscore{}.pdf'.format(savedir,cond,zscore))
            
            
        if cond =='WT': 
            
            SORTED_COEFFS = np.asarray(R_coeffs)
            SORTED_P = np.asarray(P_values)
            
        else : 
            SORTED_COEFFS = np.vstack((SORTED_COEFFS, np.asarray(R_coeffs)))
            SORTED_P = np.vstack((SORTED_P, np.asarray(P_values)))
    
    except ZeroDivisionError :
        continue

Rdf = pd.DataFrame(SORTED_COEFFS, index=conditions, columns=R_labels)
PdF = pd.DataFrame(SORTED_P, index=conditions, columns=R_labels)

if savedata == True : 
    
    Rdf.to_excel('{}/00_SORTED_COEFFICIENTS_MICROZONES.xlsx'.format(savedir))
    PdF.to_excel('{}/00_SORTED_PVALUES_MICROZONES.xlsx'.format(savedir))

    