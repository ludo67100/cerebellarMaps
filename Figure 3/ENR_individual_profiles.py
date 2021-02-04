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
import math
from scipy import stats


file = 'D:/000_PAPER/Enrichissement/Enrichment.xlsx'


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
    
    print(a,b)

    # compute r^2
    r2 = 1.0-(sum((y-f(x,a,b))**2)/((n-1.0)*np.var(y,ddof=1)))
    
    # calculate parameter confidence interval
    a,b = unc.correlated_values(popt, pcov)

    
    # plot data
    if plot == True:
        plt.figure()
        plt.scatter(x, y, s=20, label='Data',color='0.5',alpha=0.3)
    
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


if __name__ == '__main__':

    
    
    dataset = pd.read_excel(file,header=0,index_col=0,na_rep='nan')
    
    plt.figure()
    plt.title('ENR individual plots')
    
    plt.ylabel('Distance per hour (m)')
    
    profiles = []
    
    X, Y = [],[]
    
    for animal in dataset.index: 
        
        profile = dataset.loc[animal,'Day1':'Day19']
        
        #print(profile)
        
        plt.plot(profile, label=animal,color='0.5',alpha=0.3)
        
        profiles.append(profile.values)
        
        
        for i in range(len(profile)):
            
            if math.isnan(profile[i]) == False:
                
                Y.append(profile[i])
                X.append(i+1)
    
    
    Avg = np.nanmean(profiles,axis=0)
    Std = np.nanstd(profiles, axis=0)
    
    plt.plot(Avg, linewidth=2,color='green')
    plt.fill_between(np.arange(0,len(Avg),1),Avg+Std, Avg-Std, alpha=0.2,color='green')
    
    LinReg(np.asarray(X),np.asarray(Y))
    
    print (stats.linregress(X,Y))
    
        
    #plt.legend(loc='best')

    

