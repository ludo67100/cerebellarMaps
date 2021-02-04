# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:52:52 2020

@author: Ludovic.spaeth
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
import math 
import scipy.stats as sp

def deal_with_nans(x,y):
    '''
    Removes nans in 1D pattern
    returs truncated pattern and corresponding position array
    '''
    
    assert len(x) == len(y), 'x & y do not have the same length'
    
    new_x, new_y = [],[]
    for i in range(len(y)):
        
        if math.isnan(y[i]) == False:
            new_x.append(x[i])
            new_y.append(y[i])
            
    return np.asarray(new_x), np.asarray(new_y)
        

def Create_Continuous_Y_Array(x,y):
    '''
    Interpolate discrete pattern 
    '''
    #TODO : This scaling assume a step of one µm/unit per point.
    convolvedx=range(int(min(x)),int(max(x))+1)
    interpolatedy=np.zeros(len(convolvedx))
    #TODO : This transform a point to its nearest integer. So it introduces an error with small numbers
    x=[int(i) for i in x]
    #set a value at 
    counter=0
    for j,i in enumerate(convolvedx):
        if int(i) in x:
            interpolatedy[j]=y[counter]
            counter+=1 
     
    return np.array(convolvedx),np.array(interpolatedy)   
        

def Convolve(y,
             x=None,
             func='gaussian',
             sigma=1.,
             mean=0.,
             Range=(None,None),
             Normalize=False,
             Rescaled=False,
             Plateaued=False,
             ConvolvedTraceRendering=False,
             ShowFunc=False):
    '''
    -The function convolve the discrete y array with a gaussian array of 'sigma' sd.
    -an x axis can be defined (must be the same length than y)
    -Range is a tuple cropping the data between minx and maxx (or any defined limits)
     The Range parameter prevent border effects due to the convolution
    -Normalize set the biggest peak at 1
    -Rescale readjust the peak size to fit original data (override Normalize = True) 
    -Plateaued = True is necessary when using binnarized mode, else Plateaued should be False
    '''
    if Range==(None,None):
        Min=int(min(x))
        Max=int(max(x))
    else:
        Min=Range[0]
        Max=Range[1]
        

    ############# Values are prepared for convolution ###############
    #The axis in which all line average will be done. It has to be the same for all averages
    #so we use global values for that 
#    global LowerMapLimit
#    global HigherMapLimit
#    global ShowConvolutionModel
    
    axis = np.linspace(Min,Max,Max-Min) 
    
    Interpolationfunction=ReturnFunction(axis,mean=mean,sigma=sigma,func=func)

#    if ShowFunc == True and ShowConvolutionModel==True:
#        RenderConvolutionFunction(axis,Interpolationfunction,func)
        
    convolvedx,interpolatedy=Create_Continuous_Y_Array(x,y)    
    
    minIdx = np.where(convolvedx==Min)[0][0]
    maxIdx = np.where(convolvedx==Max)[0][0]
    
    
    #############Convolution###############
    if len(Interpolationfunction)<len(interpolatedy):
        convolvedy=np.convolve(interpolatedy, Interpolationfunction, mode='same')
    else:
        # TODO : when diff is shorter that Interpolationfunction, there is a probleme here.
        #I could only fixit by shortening the convolution function. Ideal interpolatedy should be extended with NaNs
        diff=int(len(Interpolationfunction)-len(interpolatedy))
        Interpolationfunction=Interpolationfunction[diff//2:-diff//2]
        convolvedy=np.convolve(interpolatedy, Interpolationfunction, mode='same')
        
    return convolvedx[minIdx:maxIdx],convolvedy[minIdx:maxIdx]

def ReturnFunction(axis,func='gaussian',mean=0.,sigma=1.):
    '''
    This function create the convolution function
    Convolution function will be of the same length than axis
    '''
    
    if func == 'gaussian':
        Curve=sp.norm.pdf(axis, mean, sigma)
        Curve/=np.max(Curve)
    elif func == 'tophat':
        Curve=np.zeros(len(axis))
        Curve[len(axis)/2-int(sigma):len(axis)/2+int(sigma)]=1
    elif func == 'cosine':
        Curve=2*np.cos(axis/sigma)
        Curve[0:len(axis)/2-int(sigma*np.pi/2.)]=0
        Curve[len(axis)/2+int(sigma*np.pi/2.):]=0
        Curve/=np.max(Curve)
    elif func == 'Rcosine':
        #TODO: Check if real RCosine
        Curve=np.cos(axis/sigma)
        Curve/=np.max(Curve)*0.5
        Curve=np.array([n if n <1. else 1. for n in Curve])
        Curve[0:len(axis)/2.-int(sigma*np.pi/2.)]=0.
        Curve[len(axis)/2.+int(sigma*np.pi/2.):]=0.
    elif func == 'triangle':
        centre_width = int(2*sigma)
        side_width = int((len(axis) - centre_width)/2)
        sigma = int(sigma)
        Curve = np.hstack((np.zeros(side_width), np.linspace(0,sigma,sigma), np.linspace(sigma,0,sigma), np.zeros(side_width)))       
        Curve/=max(Curve)
    else:
        print("Function parameters not recognised, Script aborted")
        exit() 
    
    return Curve


def column_activation(zscore, cutoff=1.96):
    '''
    Computes column activation on a 2D map
    Returns list of activation (in %) along the mediolateral axis
    '''
    assert type(zscore) == np.ndarray, 'Input Zscore is not an array'
    
    profile = []
    for i in range(zscore.shape[1]):
        
        count = round(float(len([x for x in zscore[:,i] if x >= cutoff])) / float(zscore.shape[0]), 3)
        
        profile.append(count)
        
    return profile


def range_patterns(x,y,xrange=[-200,200]):
    '''
    Refile or truncate patterns to fit them in the range
    Prior step for mediolateral alignement and/or median analysis
    '''
    xstep = x[1]-x[0]
    newx = np.arange(xrange[0],xrange[1]+xstep,xstep)
       
    newy=[]

        
    for i in range(len(newx)):
        
        position = newx[i]
        
        if position in x:
            newy.append(y[np.where(x==position)[0][0]])
        else:
            newy.append(np.nan)
    
    return newx, newy
    

def SEM(data,axis=0):
    '''
    Returns Standard Error to the Mean for a given distribution
    '''
    return np.nanstd(data,axis=axis)/np.sqrt(len(data))

def MAD(a,axis=None):
    '''
    Computes median absolute deviation of an array along given axis
    '''
    #Median along given axis but keep reduced axis so that result can still broadcast along a 
    
    med = np.nanmedian(a, axis=axis, keepdims=True)
    mad = np.nanmedian(np.abs(a-med),axis=axis) #MAD along the given axis
    
    return mad 


def bootstrap_patterns(patterns, run=10000, N=10, input_method='median', output_method='median'):
    
    '''
    Bootstraps synaptic patterns and returns median or average pattern
    
    patterns (list of arrays) : the patterns (data)
    run (int) : the amount of runs for average/median
    N (int) : number of draws for each cycle
    input_method (str) : 'median' , 'average' stores median or average value for each run 
    output_method (str) : 'median' or 'average' : returns medianed or averaged pattern
    
    '''
    
    endCycle = []
    
    for i in range(run): 
        
        temp = []
        
        for j in range(N):
        
            randIndex = np.random.randint(0,len(patterns),size=1)[0]
                        
            temp.append(patterns[randIndex])
            
            if len(temp) == N:
                pass
            else:
                continue
                
        if input_method == 'median' : 
            endCycle.append(np.nanmedian(temp, axis=0))
        
        elif input_method == 'average' : 
            endCycle.append(np.nanmean(temp, axis=0))


    if output_method == 'median': 
        out_bootstrap = np.nanmedian(endCycle, axis=0) 
        out_deviation = MAD(endCycle, axis=0)
        
    elif output_method == 'average': 
        out_bootstrap = np.nanmean(endCycle, axis=0)
        out_deviation = np.nanstd(endCycle, axis=0)
        
    return np.asarray(out_bootstrap), np.asarray(out_deviation)


def smooth_curve(y, window_size=51, polynomial_order=3):
    '''
    Returns smoothed y curve with savgol filter
    y (1D-array) : the data
    window_size (int) : The length of the filter window (i.e. the number of coefficients). window_length must be a positive odd integer
    polynomial_order (int) : he order of the polynomial used to fit the samples. polyorder must be less than window_length
    '''
    
    from scipy.signal import savgol_filter
    
    return savgol_filter(y, window_size, polynomial_order)
    

#-----------------------------------------------------------------------------------------------------------------    

#Define the groups
groups = ['P9P10','P12P13','P14P18','P30P40']
colors = ['lightskyblue','skyblue','deepskyblue','royalblue']




#Input folder - AMPLITUDE-----------------------------------------------------
dataSource = 'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS'

#Savedir 
saveDir = 'E:/000_PAPER/Development/Amplitude_Analysis/01_MEDIAN'

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'
positionFileType = 'Positions_cp_centered_OK.csv' 

#
###Input folder - CHARGE---------------------------------------------------------
#dataSource = 'E:/000_PAPER/Charge_Analysis/Sigma'
#
##Savedir 
#saveDir = 'E:/000_PAPER/Charge_Analysis/Sigma/01_MEDIAN'
#
##Target file type
#fileType = 'Charge_2D_OK.csv'
#zscoreFileType = 'Charge_zscore_2D_OK.csv'
#positionFileType = 'Positions_cp_centered_OK.csv' 


zscoreCut = 3.09

#Sigma for convolution
sigma = 9

#Range in %P1- for refiling and median analysis 
mapRange = [-200,220]

#Do we save the data ?
saveData = True

#Display plot with raw, interpolated and convovled trace for each experiment (lots of plots)
showPlot = False

#Do we save the figure ?
saveFig = True

#Do we perform the stats ?
statsToDo = False


#FOR THE WHOLE MAP
MEASURES,CONDITIONS = [],[]

for group,index in zip(groups,range(len(groups))):
    print('Group = {}'.format(group))
    
    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir)]
    
    fig, ax = plt.subplots(1,2 ,figsize=(20,5))
    plt.suptitle(group)
    ax[0].set_title('Convolved patterns (sigma={})'.format(sigma))
    ax[0].set_xlabel('Mediolateral axis (%P1-)')
    ax[0].set_ylabel('Zscore')
    
    ax[1].set_title('Group Data'.format(sigma))
    ax[1].set_xlabel('Mediolateral axis (%P1-)')
    ax[1].set_ylabel('Median Zscore')    
    
    ax[1].set_ylim(0,50)
       
    allPositions, allPatterns = [],[]
    allPositionsActivation, allActivations = [],[]
    for manip,idx in zip(listOfExperiments,range(len(listOfExperiments))):
        
        print(manip)

        #Get amplitudes or charges
        measures = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,fileType),delimiter=','))
        #Get corresponding zscores
        zscores = np.abs(np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,zscoreFileType),delimiter=','))
        #Get positions 
        positions = np.genfromtxt('{}/{}/{}_{}'.format(inputDir,manip,manip,positionFileType),delimiter=',')
        
        #Column activation 
        activation = column_activation(zscores,cutoff=zscoreCut)
        
        #Amplitude & Zscore patterns
        oneDzscore = np.nanmax(zscores,axis=0)
        oneDpattern = np.nanmax(measures,axis=0)
        
        #Interpolate, convolve and align patterns
        cleand_x, cleaned_y = deal_with_nans(positions,oneDzscore)
        continuousX, continuousY = Create_Continuous_Y_Array(cleand_x,cleaned_y)
        convolveX,convolveY = Convolve(continuousY,x=continuousX,func='triangle',sigma=sigma,Range=(None,None))
        alignedx, alignedy = range_patterns(convolveX,convolveY,xrange=mapRange)
        
        ax[0].plot(alignedx, alignedy,label=manip)

        allPositions.append(alignedx)
        allPatterns.append(alignedy)

        if showPlot == True:
            #Do a figure
            plt.figure()
            plt.title('{} Zscore Pattern'.format(manip))
            plt.ylabel('Zscore')
            plt.xlabel('Position on the mediolateral axis')
            plt.plot(cleand_x,cleaned_y,label='Raw',linestyle='solid')
            plt.plot(continuousX,continuousY,label='Interpolated',linestyle='dashed')
            plt.plot(convolveX,convolveY,label='Convolved (sigma={} microns)'.format(sigma),linestyle='solid')
            plt.plot(alignedx,alignedy,label='Aligned'.format(sigma),linestyle='dashed')
            plt.legend(loc='best')
            
    
    xaxis = np.arange(mapRange[0],mapRange[1]+1,1)        
    
    
    #Median Zscore
    medianZscore = np.nanmedian(allPatterns,axis=0)
    MaD = MAD(allPatterns,axis=0)

    ax[1].fill_between(xaxis,0,medianZscore+MaD,color='0.8',label='MAD')
    ax[1].fill_between(xaxis,0,medianZscore,color='0.2',label='ns')
    ax[1].fill_between(xaxis,zscoreCut,medianZscore,where=medianZscore>=zscoreCut,color=colors[index],interpolate=True,label='connected')
    
    #Bootstrap
    zeBootsrap, deviation = bootstrap_patterns(allPatterns,run=10000, N=len(allPositions),input_method='median', output_method='average')
    
    #zeBootsrap = smooth_curve(zeBootsrap, 21,3)
    
    #Plot bootstrap curve way
    ax[1].plot(xaxis,zeBootsrap, label='Bootstrap', color='red')
    ax[1].fill_between(xaxis, zeBootsrap-deviation, zeBootsrap+deviation, alpha=0.2,color='red',label='SD')
    
    #Zscore limit
    zscoreLine = np.ones(len(xaxis))*zscoreCut
    ax[1].plot(xaxis, zscoreLine, linestyle='--', color='black', label='Threshold')
    
     #Plot bootstrap filled way
#    ax[1].fill_between(xaxis,0,zeBootsrap+deviation,color='0.8',label='SD')
#    ax[1].fill_between(xaxis,0,zeBootsrap,color='0.2',label='ns')
#    ax[1].fill_between(xaxis,zscoreCut,zeBootsrap,where=zeBootsrap>=zscoreCut,color=colors[index],interpolate=True,label='connected')
    

    #Add the zebrin bands
    for i in range(2):
        ax[i].axvspan(xmin=-140.5,xmax=-114.6,color='green',alpha=0.2)
        ax[i].axvspan(xmin=-11.5,xmax=0,color='green',alpha=0.2)
        ax[i].axvspan(xmin=100,xmax=125.5,color='green',alpha=0.2)
        ax[i].legend(loc='best')
    
    
    if saveFig == True: 
        fig.savefig('{}/{}_MedianPattern.png'.format(saveDir, group))
        fig.savefig('{}/{}_MedianPattern.pdf'.format(saveDir, group))


#Check random draw
#plt.figure()
#draws=[]
#for i in range(10000):
#    for j in range(14):
#        randIndex = np.random.randint(0,14,size=1)[0]
#        draws.append(randIndex)
#        
#plt.hist(draws, bins=100, width=0.8)
#plt.xlabel('Pattern Index (draw)')


        
        
        
        
     
