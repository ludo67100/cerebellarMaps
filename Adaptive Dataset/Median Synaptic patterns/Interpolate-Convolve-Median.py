# -*- coding: utf-8 -*-
"""
Created on Sun May 24 13:58:44 2020

Computes median synaptic patterns (high resolution) for each behavioral condition

@author: ludov
"""

#Define the groups
groups = ['WT','ENR1','ENR2','EC','ES','LC','LS']
newNomConditions = ['control','short training','long training','early cuff','early sham','adapted cuff','adapted sham']
colors = ['skyblue','lightgreen','green','orange','purple','lightcoral','0.5']

#Input folder - AMPLITUDE-----------------------------------------------------
#dataSource = 'C:/Users/ludov/Documents/Spaeth_Bahuguna_et_al/ANALYSIS/DATASET/EPHYS/Adaptive_Dataset'
dataSource = 'D:/03_FORMATED_DATA/For Paper/EPHYS/Adaptive_Dataset'
#dataSource = './For Paper/EPHYS/Adaptive_Dataset'

#Savedir 
saveDir = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/REVISION CODE/CODE/FIGURE 4/C'
#saveDir = './data/'


zscoreCut = 3.09

#Sigma for convolution
sigma = 9

#Range in %P1- for refiling and median analysis 
# mapRange = [-200,200]
mapRange = [-200,200]

#Do we save the data ?
saveData = False

#Display plot with raw, interpolated and convovled trace for each experiment (lots of plots)
showPlot = True

#Do we save the figure ?
saveFig = False

#Do we compile the individual patterns for later RF analysis ?
compilePatternsForRF = False

#Do we do boostrap on patterns ? Takes a few seconds more to compute
do_bootstrap = False

#Do we perform the stats ?
statsToDo = True

#Method for cumulative : 'median' or 'mean' 
cumMethod = 'mean'

#For average, use SD or SEM 
averageSDorSEM = 'SD'

#Replace non significant columns with 0 for cumulatives 
replaceSilentColumnsWithZeros = False

#Store the MADs for comparison 
storeMad = []

#Kwards for KS test
ksMode = 'auto'
ksAlt = 'two-sided'

#Method for KS p values correction
multiCorrMethod = 'holm'



import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import numpy as np
import matplotlib.pyplot as plt
import os 
import math 
import scipy.stats as sp
import pandas as pd 
from scipy import stats
import warnings
from statsmodels.stats.multitest import multipletests as multi

#Ignore nan-slice warnings
warnings. filterwarnings(action='ignore')


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
        
        count = round(float(len([x for x in zscore[:,i] if x >= cutoff])) / float(zscore.shape[0]), 3)*100
        
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
    return stats.sem(data, axis=axis, nan_policy='omit')

def MAD(a,axis=None):
    '''
    Computes median absolute deviation of an array along given axis
    '''
    #Median along given axis but keep reduced axis so that result can still broadcast along a 
    
    # med = np.nanmedian(a, axis=axis, keepdims=True)
    # mad = np.nanmedian(np.abs(a-med),axis=axis) #MAD along the given axis
    
    # return mad 
    return stats.median_abs_deviation(a,nan_policy='omit',axis=0)


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


def bootstrap_patterns_random(patterns, run=10000, N=10, input_method='median', output_method='median'):
    
    '''
    Bootstraps synaptic patterns and returns median or average pattern. The pattern pot is randomized at each turn 
    
    patterns (list of arrays) : the patterns (data)
    run (int) : the amount of runs for average/median
    N (int) : number of draws for each cycle
    input_method (str) : 'median' , 'average' stores median or average value for each run 
    output_method (str) : 'median' or 'average' : returns medianed or averaged pattern
    
    '''
    
    endCycle = []
    
    for i in range(run): 
        
        patterns = shuffle_these_patterns(patterns)
        
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


def shuffle_these_patterns(listOfPatterns):
    '''
    Shuffle patterns (intra-patterns) from a list of patterns
    '''
    
    import sklearn
    
    listOfPatternsCopy = listOfPatterns
    
    shuffled_arrays = [sklearn.utils.shuffle(x) for x in listOfPatternsCopy]
        
    return shuffled_arrays
    

#-----------------------------------------------------------------------------------------------------------------    

#Target file type
fileType = 'Amp_2D_OK.csv'
zscoreFileType = 'Amp_zscore_2D_OK.csv'
positionFileType = 'Positions_cp_centered_OK.csv'

print('------------------------------------------------------')
print('-----------Map Range(%P1-)={}------------'.format(mapRange))
print('------------------------------------------------------') 


#FOR THE WHOLE MAP
MEASURES,CONDITIONS = [],[]

#To store cumulative patterns 
cumulativePatternsGlob, mapGlob, condGlob,PatternsGlob = [],[],[],[]

#Lists for convolved patterns and zscores
convolved_patterns, convolved_mads = [],[]
convolved_median_amps, convolved_amps_mads = [],[]
convolved_median_activation, convolved_activation_mads = [],[]

#Lists for raw amplitude profiles - will be saved in a dedicated DF for later RF 
rawAmpProfiles = []
groupForAmpProfile = []
nameForAmpProfile = []

#List for cumulative distribution at x=0 and x=max
x0Cumulatives, xMaxCumulatives = [],[]
x0CumulativesActivation, xMaxCumulativesActivation = [],[]

for group,index in zip(groups,range(len(groups))):
    print('Group = {}'.format(group))
    

    #Get input directory
    inputDir = '{}/{}'.format(dataSource,group)
    
    #Get list of experiments
    listOfExperiments = [x for x in os.listdir(inputDir) if group in x]
    
    #For ENR1 and 2
    listOfExperiments = [x for x in os.listdir(inputDir)]
    
    #Store interquartiles range for each group 
    interquartiles = []

    
    fig, ax = plt.subplots(1,2 ,figsize=(20,5))
    plt.suptitle(group)
    ax[0].set_title('Convolved patterns (sigma={})'.format(sigma))
    ax[0].set_xlabel('Mediolateral axis (%P1-)')
    ax[0].set_ylabel('Zscore')
    
    ax[1].set_title('Group Data (sigma={})'.format(sigma))
    ax[1].set_xlabel('Mediolateral axis (%P1-)')
    ax[1].set_ylabel('Median Zscore')    
    
    ax[1].set_ylim(0,25)
       
    allPositions, allPatterns, allAmplitudes = [],[],[]
    allPositionsActivation, allActivations = [],[]
    for manip,idx in zip(listOfExperiments,range(len(listOfExperiments))):
        
        #print(manip)
        mapGlob.append('{}_{}'.format(manip,newNomConditions[index])) 

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
        
        #Interpolate, convolve and align amplitudes
        cleand_xAmp, cleaned_yAmp = deal_with_nans(positions,oneDpattern)
        continuousXAmp, continuousYAmp = Create_Continuous_Y_Array(cleand_xAmp, cleaned_yAmp)
        convolveXAmp,convolveYAmp = Convolve(continuousYAmp, x=continuousXAmp,func='triangle',sigma=sigma,Range=(None,None))
        alignedxAmp, alignedyAmp = range_patterns(convolveXAmp,convolveYAmp,xrange=mapRange)
        
        #Interpolate, convolve and align column activation - Do not convolve
        cleand_xActivation, cleaned_yActivation = deal_with_nans(positions,np.array(activation))
        continuousXActivation, continuousYActivation = Create_Continuous_Y_Array(cleand_xActivation, cleaned_yActivation)
        #convolveXActivation,convolveYActivation = Convolve(continuousYActivation, x=continuousXActivation,func='gaussian',sigma=sigma,Range=(None,None))
        alignedxAactivation, alignedyAactivation = range_patterns(continuousXActivation,continuousYActivation,xrange=mapRange)

        
        if replaceSilentColumnsWithZeros==True: 
            for i in range(len(alignedy)):
                if alignedy[i] >= zscoreCut:
                    alignedyAmp[i] = alignedyAmp[i]
                else:
                    alignedyAmp[i] = 0
        else:
            pass 
        
        ax[0].plot(alignedx, alignedy,label=manip)
            

        allPositions.append(alignedx)
        allPatterns.append(alignedy)
        allAmplitudes.append(np.nancumsum(alignedyAmp))
        allActivations.append(np.nancumsum(alignedyAactivation))
        
        #Store the cumulative amplitudes and zscore patterns 
        cumulativePatternsGlob.append(alignedyAmp)
        PatternsGlob.append(alignedy)
        
        
        #Append amps in Df for RF and tSNE if maps are complete
        if math.isnan(np.sum(oneDpattern))==False:
            rawAmpProfiles.append(oneDpattern)
            groupForAmpProfile.append(group)
            nameForAmpProfile.append(manip)
            

        if showPlot == True:
            #Do a figure
            plt.figure()
            plt.title('{} Zscore Pattern'.format(manip))
            plt.ylabel('Zscore')
            plt.xlabel('Position on the mediolateral axis')
            plt.plot(cleand_x,cleaned_y,label='Raw',linestyle='solid')
            plt.plot(continuousX,continuousY,label='Interpolated',linestyle='dashed')
            plt.plot(convolveX,convolveY,label='Convolved (sigma={} microns)'.format(sigma),linestyle='solid')
            plt.plot(alignedx,alignedy,label='Aligned (sigma={})'.format(sigma),linestyle='dashed')
            plt.legend(loc='best')
            
    
    xaxis = np.arange(mapRange[0],mapRange[1]+1,1)      
    
    #Export DataFramr with aligned, convolved patterns for further analysis 
    if saveData == True:
        extDf = pd.DataFrame(allPatterns,columns=allPositions[0])
        extDf.to_excel('{}/{}_aligned_zscore_patterns.xlsx'.format(saveDir,group))
    else:
        print('saveData == False, Df was not saved')
    
    #Median Zscore
    # medianZscore = np.nanmedian(allPatterns,axis=0)
    # MaD = MAD(allPatterns,axis=0)
    
    medianZscore = np.nanmean(allPatterns,axis=0)
    MaD = np.nanstd(allPatterns,axis=0)
    
    #Get cumulative distributions x0 and xMax
    #x0Cumulatives.append(np.array(allAmplitudes)[:,int(len(medianZscore)/2.)]/1000.)
    xMaxCumulatives.append(np.array(allAmplitudes)[:,-1]/1000.)
    
    x0CumulativesActivation.append(np.array(allActivations)[:,int(len(medianZscore)/2.)])
    xMaxCumulativesActivation.append(np.array(allActivations)[:,-1])
    
    if cumMethod == 'median':
        medianAmplitude = np.nanmedian(allAmplitudes, axis=0)
        MaDamps = MAD(allAmplitudes, axis=0)
        storeMad.append(MaDamps)
        medianActivation = np.nanmedian(allActivations,axis=0)
        MaDactivation = MAD(allActivations,axis=0)
        errorMethod = 'MAD'
    else:
        medianAmplitude = np.nanmean(allAmplitudes, axis=0)

        medianActivation = np.nanmean(allActivations,axis=0)
        
        if averageSDorSEM == 'SEM':
            MaDamps = SEM(allAmplitudes, axis=0)
            storeMad.append(MaDamps)
            MaDactivation = SEM(allActivations,axis=0)
            errorMethod = 'SEM'
        else:
            MaDamps = np.nanstd(allAmplitudes, axis=0)
            storeMad.append(MaDamps)
            MaDactivation = np.nanstd(allActivations,axis=0)
            errorMethod = 'SD'

    ax[1].fill_between(xaxis,0,medianZscore+MaD,color='0.8',label='MAD')
    ax[1].fill_between(xaxis,0,medianZscore,color='0.2',label='ns')
    ax[1].fill_between(xaxis,zscoreCut,medianZscore,where=medianZscore>=zscoreCut,color=colors[index],interpolate=True,label='connected')
    
    convolved_patterns.append(medianZscore)
    convolved_mads.append(MaD)
    
    convolved_median_amps.append(medianAmplitude/1000.) #in nA instead of pA - easier for cumulative
       
    convolved_amps_mads.append(MaDamps/1000.) #in nA instead of pA - easier for cumulative

    convolved_median_activation.append(medianActivation)
    convolved_activation_mads.append(MaDactivation)
    
    if do_bootstrap == True:
    
        #Bootstrap
        zeBootsrap, deviation = bootstrap_patterns(allPatterns,run=10000, N=len(allPositions),input_method='average', output_method='average')
        
        zeBootsrap = smooth_curve(zeBootsrap, 21,3)
        
        #shuffledPatterns = shuffle_these_patterns(allPatterns)
        
        shuffleBootstrap, shuffleDeviation = bootstrap_patterns_random(allPatterns,
                                                                       run=10000, N=len(allPositions),input_method='average', output_method='average')
        
        shuffleBootstrap = smooth_curve(shuffleBootstrap,21,3)
        
        #Plot bootstrap curve way
        ax[1].plot(xaxis,zeBootsrap, label='Bootstrap', color='red')
        ax[1].fill_between(xaxis, zeBootsrap-deviation, zeBootsrap+deviation, alpha=0.2,color='red',label='SD')
        
        ax[1].plot(xaxis,shuffleBootstrap, label='Shuffle', color='blue')
        ax[1].fill_between(xaxis, shuffleBootstrap-shuffleDeviation, shuffleBootstrap+shuffleDeviation, alpha=0.2,color='blue',label='Shuffle SD')
        
    
    #Zscore limit
    zscoreLine = np.ones(len(xaxis))*zscoreCut
    ax[1].plot(xaxis, zscoreLine, linestyle='--', color='black', label='Threshold')

    #Add the zebrin bands
    for i in range(2):
        ax[i].axvspan(xmin=-140.5,xmax=-114.6,color='green',alpha=0.2)
        ax[i].axvspan(xmin=-11.5,xmax=0,color='green',alpha=0.2)
        ax[i].axvspan(xmin=100,xmax=125.5,color='green',alpha=0.2)
        ax[i].legend(loc='best')
    
    
    if saveFig == True: 
        fig.savefig('{}/{}_MedianPattern.png'.format(saveDir, group))
        fig.savefig('{}/{}_MedianPattern.pdf'.format(saveDir, group))


    if do_bootstrap == True:

        #Check random draw
        plt.figure()
        draws=[]
        for i in range(10000):
            for j in range(len(allPatterns)):
                randIndex = np.random.randint(0,len(allPatterns),size=1)[0]
                draws.append(randIndex)
                
        plt.hist(draws, bins=100, width=0.8,color=colors[index])
        plt.xlabel('Pattern Index (#)')
        plt.ylabel('Amount of draws')
        if saveFig==True:
            plt.savefig('{}/{}_PatternDraw.pdf'.format(saveDir, group))

#Cumulative amplitudes_________________________________________________________
if statsToDo == True:
    print('')
    print('--------{} Amplitude Patterns Comparison---------'.format(cumMethod))
    print('KS 2samp test - mode={} - alternative={}'.format(ksMode,ksAlt))

    #Figure for the cumulatives
    cumFig, cumAx = plt.subplots(1,len(groups),sharex=True,sharey=True,figsize=(16,4))
    cumAx[0].set_ylabel('{} cumulative Amplitude (nA+/-{})'.format(cumMethod,errorMethod))
    
    #Figure for cumulative comparison
    cumCompareFig, cumCompareAx = plt.subplots(1,1)
    cumCompareAx.set_ylabel('Cumulative[i]/Cumulative[CTRL]')
    cumCompareAx.set_xlabel('Distance (%P1-)')

    
        
    risePlot, riseAx = plt.subplots(1,1,sharex=True,sharey=True)
    riseAx.set_ylabel('cumulative Amplitude (nA +/-{})'.format(errorMethod))

    
    assert groups[0] == 'WT', 'the ref pattern for ks test is not the CTRL pattern'
    ref_median_pattern = convolved_median_amps[0]
    
    ks_pvalues = []
    for med, mad, group, i in zip(convolved_median_amps,convolved_amps_mads,groups,range(len(groups))):
        
        ks_test = sp.ks_2samp(ref_median_pattern,med,mode=ksMode,alternative=ksAlt)
        print('    KS test WT vs {} stat = {}'.format(group,ks_test[0]))
        print('    KS test WT vs {} p val = {}'.format(group, ks_test[1]))
                
        ks_pvalues.append(ks_test[1])

        cumAx[i].plot(xaxis,convolved_median_amps[0],color=colors[0],label=groups[0])
        cumAx[i].fill_between(xaxis,convolved_median_amps[0]+convolved_amps_mads[0],convolved_median_amps[0]-convolved_amps_mads[0],color=colors[0],alpha=.2)
        
        cumCompareAx.plot(xaxis,convolved_median_amps[i]/convolved_median_amps[0],color=colors[i],label=groups[i])
        cumCompareAx.legend(loc='best')
        
        cumAx[i].plot(xaxis,med,color=colors[i],label=groups[i])
        cumAx[i].fill_between(xaxis,med+mad,med-mad,color=colors[i],alpha=.2)

  
        cumAx[i].set_xlabel('Distance (%P1-)')
        cumAx[i].legend(loc='best')
        
        maxRiseAmp = med[-1] 
        maxRiseMad = mad[-1]
        
        #Get values at x=0
        # x0 = [x for x, y in zip(range(len(xaxis)), xaxis) if y == 0][0]
        # midRiseAmp = med[x0]
        # midRiseMad = mad[x0]

        #print('        at x=0. {} Cumulative Amplitude in {} (nA+/-{})= {}+/-{}'.format(cumMethod,group,errorMethod,round(midRiseAmp,2),round(midRiseMad,2)))
        
        print('        Max. {} Cumulative Amplitude in {} (nA+/-{})= {}+/-{}'.format(cumMethod,group,errorMethod,round(maxRiseAmp,2),round(maxRiseMad,2)))
        print('')
        #box50 = riseAx[0].boxplot(x0Cumulatives,vert=False,patch_artist=True,showmeans=True)
        box100 = riseAx.boxplot(xMaxCumulatives,vert=True,patch_artist=True,showmeans=True)
    
    #Apply multitest corrections on KS p values
    checkKSpval = multi(ks_pvalues,alpha=0.05, method=multiCorrMethod, is_sorted=False)
    print('Correction for multiple KS comparison')
    print(checkKSpval)
    
    for p in range(len(checkKSpval[1])):
        print ('Corrected ({}) p-value for WT vs {} (KS)= {:.4g}'.format(multiCorrMethod,groups[p],checkKSpval[1][p]))
        cumAx[p].set_title('p={:.3g}'.format(checkKSpval[1][p]))   
    
    #Apply colors to boxplots
    for patch100, color in zip(box100['boxes'],colors):
        patch100.set_facecolor(color)
        
        
    #Test xMaxCumulative distributions--------------------------------------------
    kruskalMax = stats.kruskal(xMaxCumulatives[0],
                               xMaxCumulatives[1],
                               xMaxCumulatives[2],
                               xMaxCumulatives[3],
                               xMaxCumulatives[4],
                               xMaxCumulatives[5])
    
    print('--- Comparison on x=Max ---') 
    print(kruskalMax)
    riseAx.set_title('Distributions at x=max - p(KW)={:.4g}'.format(kruskalMax[1]))

    print ('    Post Hoc MWU/T test on CTRL vs conditions')
    for k in range(len(xMaxCumulatives)):
        
        normality = stats.shapiro(xMaxCumulatives[k])
        if normality[1] < 0.05: 
            print('{} distribution does not fit normality'.format(groups[k]))
        else:
            print ('{} distribution follows normality'.format(groups[k]))
            
        if stats.shapiro(xMaxCumulatives[0])[1] >= 0.05:
            
            if stats.shapiro(xMaxCumulatives[k])[1] >=0.05:
                postHocTtest = stats.ttest_ind(xMaxCumulatives[0],xMaxCumulatives[k])
                print ('    {} group vs {} group (t test)'.format(groups[0],groups[k]))
                print ('    {}'.format(postHocTtest))
                print ('')
                
            else:
                postHocMwu = stats.mannwhitneyu(xMaxCumulatives[0],xMaxCumulatives[k])
                
                print ('    {} group vs {} group (MWU)'.format(groups[0],groups[k]))
                print ('    {}'.format(postHocMwu))
                print ('')
        else:
            postHocMwu = stats.mannwhitneyu(xMaxCumulatives[0],xMaxCumulatives[k])
            
            print ('    {} group vs {} group (MWU)'.format(groups[0],groups[k]))
            print ('    {}'.format(postHocMwu))
            print ('')
            

if compilePatternsForRF == True:
    #-------------Compile DataFrames for later RF analysis------------
    #One DF with sorted labels, then another with shuffled labels 
    rawAmpProfilesDf = pd.DataFrame(rawAmpProfiles)
    rawAmpProfilesDf['Condition'] = groupForAmpProfile
    rawAmpProfilesDf['name'] = nameForAmpProfile
    rawAmpProfilesDf.to_excel('{}/RawAmpProfilesForRandomForest_SORTED_LABELS.xlsx'.format(saveDir))
            
            
    import random 
    #Fix the see
    random.seed(4)
    rawAmpProfilesDfShuffle = pd.DataFrame(rawAmpProfiles)
    rawAmpProfilesDfShuffle['Condition'] = random.sample(groupForAmpProfile,len(groupForAmpProfile))
    rawAmpProfilesDfShuffle['name'] = random.sample(nameForAmpProfile,len(nameForAmpProfile))
    rawAmpProfilesDfShuffle.to_excel('{}/RawAmpProfilesForRandomForest_RANDOM_LABELS.xlsx'.format(saveDir))
    
    
#Save the dataframes -----------------------------------------------------------

cumulativeGlobDf = pd.DataFrame(data=cumulativePatternsGlob,index=mapGlob).T
cumulativeGlobDf.to_excel('{}/Cumulative_Amplitude_Patterns.xlsx'.format(saveDir))
patternsDf = pd.DataFrame(data=PatternsGlob, index=mapGlob,columns=xaxis).T
patternsDf.to_excel('{}/Aligned_Zscore_Patterns.xlsx'.format(saveDir))



#------------------Compared each condition early vs adapted ------------------------
groups = ['WT','ENR1','ENR2','EC','ES','LC','LS']

#Early vs adapted Cuff
earlyCuffCum = convolved_median_amps[3]
adaptedCuffCum = convolved_median_amps[5]

adaptedVsEarlyCuff = sp.ks_2samp(earlyCuffCum,adaptedCuffCum,mode=ksMode,alternative=ksAlt)
print('Early cuff vs adapted cuff comparison (KS, mode={}, alt={}) = {}'.format(ksMode, ksAlt,adaptedVsEarlyCuff))
print()

#Early vs adapted Sham
earlyShamCum = convolved_median_amps[4]
adaptedShamCum = convolved_median_amps[6]

adaptedVsEarlySham = sp.ks_2samp(earlyShamCum,adaptedShamCum,mode=ksMode,alternative=ksAlt)
print('Early sham vs adapted sham comparison (KS, mode={}, alt={}) = {}'.format(ksMode, ksAlt,adaptedVsEarlySham))
print()

#short vs long train
shortTrainCum = convolved_median_amps[1]
longTrainCum = convolved_median_amps[2]

shortVsLongTraining = sp.ks_2samp(shortTrainCum,longTrainCum,mode=ksMode,alternative=ksAlt)
print('short vs long training comparison (KS, mode={}, alt={}) = {}'.format(ksMode, ksAlt,shortVsLongTraining))
print()


        
     
