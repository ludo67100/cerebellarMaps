# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 18:24:00 2020

@author: Ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 17:46:35 2019
---------------AMPLITUDE WAY------------------------------------------------------------------
Despite all my efforts, this scripts takes a wile (~2-3 minutes) to compute and build a full map. Sorry. 
UPDATE 06/05/2019 : This version runs with neo 0.8 (dev) and rawio loading method
please be sure to use the neo version I've modified to correct the sampling rate issue with WCP
If false, go to //site-packages/neo/rawio/winwcprawio.py, REPLACE LINE 81:
    
    assert np.unique(all_sampling_rate_interval).size == 1
    
    by
    
    self._sampling_rate = 1. / all_sampling_interval[1]
This will attribute correct sampling rate to all segments in the recording
Script built to compute Synaptic maps from in vitro uncaging experiments, in a "automatic" way
What the script does : 
    - loads WinWcp (.wcp) files containing voltage clamp recordings
    - Gets info for spacing, orientation and so on from GC_PC_mappings_data_info.xlsx excel file
    - Gets Zebrin positions and info from Mesures_ZII_LowRes_Adult_and_Dev.xlsx excel file
    
    - Calculate average amplitude, average noise, and zscore from ephy data
    - Build synaptic map in consensus orientation 
    
    - Save figures and tables  (if you want to)
@author: ludovic.spaeth
"""

from neo.rawio import WinWcpRawIO as winraw
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import pandas as pd

#-----------------------DIRECTORY & INFO---------------------------------------
#--------------CHECK ALL THESE INFO BEFORE RUNNING----------------------------- 
#------------------------------------------------------------------------------

DONE = [] #To append names of maps that are done 

NAMES = ['151216(1)','151216(2)',
         '151217(1)','151217(2)','151217(3)',
         '151128(1)',
         '151221(1)','151221(2)',
         '151117(1)',
         '151222(1)','151222(2)',
         '151229(1)','151229(2)',
         '151223(1)','151223(2)']

FOLDERS = ['151216(1)','151216(2)',
         '151217(1)','151217(2)','151217(3)',
         '151128(1)',
         '151221(1)','151221(2)',
         '151117(1)',
         '151222(1)','151222(2)',
         '151229(1)','151229(2)',
         '151223(1)','151223(2)']

NAMES = ['151117(1)']

FOLDERS = ['151117(1)']

#Where are the data ? The data are in the kitchen :) 
datadir = 'U:/RAW DATA/Development Dataset/P14P18'

#Where to save the data 
savedir = 'E:/000_PAPER/Development/Amplitude_Analysis/00_MAPS/P14P18'

#Excel sheet with zebrin band measures : Mesures_ZII_HighRes_WT_Cuff_Sham_Enr.xlsx
df  = pd.read_excel('E:/000_PAPER/Mesures_ZII_LowRes_Adult_and_Dev.xlsx',
                    sheet_name='DEV',index_col=0,header=1)  #DataFrame containing zebrin band file 

#Excel sheet with files and parameters : GC-PC_mappings_data_info.xlsx  
filelist  = pd.read_excel('E:/000_PAPER/GC-PC_mappings_data_info.xlsx',
                          sheet_name='DVLPMT',index_col=0,header=0)  #DataFrame containing files

#Number of bins for noise fitting, 100 if map is complete, 50 for security
number_bin_fit = 30

SAVING = True #If True, saves the data in savedir
CLOSE_FIG = True #If true, will close figures after computing 
FAST = False #If True, will NOT plot/save current plots to increase speed 

#For a single experiment, u can give the sigma value (if the fit fails)
#To skip and rely on fit, apply None
# Here 8.01 is applied to 151117(1)
aPrioriSigmaValue = None

#----------------------NOW THE SERIOUS BUSINESS--------------------------------
#--------------------------DO NOT MODIFY---------------------------------------
#--------------Iterate + map computation for each experiment-------------------

for NAME,FOLDER,IDX in zip(NAMES,FOLDERS,range(len(NAMES))):
    
    try : 
    
        name = NAME
        manip = NAME.replace('WT_','') #For WT nomenclature only 
        
        #Path to the data
        directory = '{}/{}'.format(datadir,FOLDER)
    
        #Channel ID     
        ch = int(filelist.loc[[manip],['Channel']].values.ravel()[0]) #Channel 0 for Im1 and 2 for Im2
    
        #Gets scanspot list order     
        reshape = int(filelist.loc[[manip],['Ordre']].values.ravel()[0])
        global_reshape = np.array([int(x) for x in str(reshape)])
    
        #-------------------------GETS FILELIST----------------------------------------
      
        print ('---------------------------------')
        print ('exp {} out of {}'.format(IDX+1,len(NAMES)+1))
        print ('Collecting files for {} experiment'.format(name))
        
        ZSOCRE_CUT = 2.0
        
        #•Get protocol type 
        protocol = str(filelist.loc[[manip],['Protocol']].values.ravel()[0])
        print (protocol)
    
        global_reshape = global_reshape-1
        
        scanspot_list = ['Scanspot 1','Scanspot 2','Scanspot 3','Scanspot 4']
        
        files = []
        
        for i in range(global_reshape.size):
            idx_list = np.ravel(filelist.loc[[manip],[scanspot_list[i]]].values)
            record_list  = []
            
            for idx in idx_list :
                if isinstance(idx,str) == True:
        
                    record_list.append(idx)
                    
            files.append(record_list)
            
        files = np.asarray(files)
        
        #----------------------ZEBRIN BANDS & ORIENTATION-------------------------------
        _BANDS_micron =  df.loc[['%s'%name],'P2- contra':'P2- ipsi']  #Selects values in xcl file from P2- contra to P2- ipsi
            
        _BANDS_norm =  df.loc[['%s norm_P1-'%name],'P2- contra':'P2- ipsi']  #Selects values in xcl file from P2- contra to P2- ipsi
        
        orientation= df.loc[['%s'%name],['Position']]   # Gets the piece of dataframe
        orientation= orientation.iloc[0,0]              # Gets the float value in the df 
        orientation= int(orientation)                   #Converts float in integer 
            
        position_left = df.loc[['%s'%name],['Pos Left']] # 32 if cell has been centered, may change in case of paired recordings 
        position_left = position_left.iloc[0,0]
        position_left = int(position_left)
            
        position_right = df.loc[['%s'%name],['Pos Right']] # 32 if cell has been centered, may change in case of paired recordings 
        position_right = position_right.iloc[0,0]
        position_right = int(position_right)
        
        
        P2mcontra = _BANDS_norm.iloc[0,0]
        P2pcontra = _BANDS_norm.iloc[0,1]
        P1mcontra = _BANDS_norm.iloc[0,2]
        P1p	= _BANDS_norm.iloc[0,3]
        BordgaucheP1m = 0.
        BorddroitP1m = 100.
        cp = _BANDS_norm.iloc[0,5]
        P2pipsi   = _BANDS_norm.iloc[0,6]
        P2mipsi   = _BANDS_norm.iloc[0,7]
                
        zebrin = np.array([P2pcontra,P1mcontra,P1p,BordgaucheP1m,cp,BorddroitP1m,P2pipsi,P2mipsi])
        
        
#        #High rest grid
#        grid = np.array([[1,17,5,21,9,25,13,29,2,18,6,22,10,26,14,30],
#                         [50,35,54,39,96,42,91,46,79,36,53,40,85,45,56,48],
#                         [31,3,19,7,23,11,27,15,32,4,20,8,24,12,28,16],
#                         [59,66,89,60,75,58,84,71,88,62,95,77,82,57,74,87],
#                         [33,73,37,70,41,65,44,80,34,68,38,93,43,94,47,61],
#                         [64,49,78,52,76,69,86,63,92,51,72,55,83,67,90,81]])
    
        #Low Res grid
        grid = np.array([[1,8,13,29,17,26,10,4],
                         [16,25,31,5,21,30,19,15],
                         [7,20,28,14,32,24,6,22,],
                         [3,11,23,9,18,27,12,2]])
                         
        sites = grid.size
        
        gridX = grid.shape[0]
        gridY = grid.shape[1]
        

        frameY = gridY*4
        
        #Position et bords de la carte ----------------------------------------------------------------------------------------------
        print ('PC location: ',cp,"%",'of P1-'	)						
        p1_minus = _BANDS_micron.iloc[0,4]			                               # Taille de P1- en microns
        print ('P1- length: ' ,p1_minus, ' microns')
        step = float(40)														# Taille de sites de photostim en microns
        print (step,' microns/photstimulation site')
        step2 = float((step/p1_minus)*100)										# Calcul du step en % de P1-
        print ('hence ',step2," %",'de P1-')
        lb = position_left*step/p1_minus*100									# Calcul bord gauche de la map en % de P1-
        print ('Map left limit = ',lb,"%")
        rb = position_right*step/p1_minus*100								    # Calcul du bord droit de la map en % de P1-
        print ('Map right limit = ',rb,"%")
        
        step_positions = (lb+rb)/(frameY-1)
        positions = np.linspace(-lb,rb,frameY)
        positions_cp_centered = np.linspace(-lb,rb,frameY)+cp
        

        
        
        _GRIDS_OF_SIGNAL_AMPLITUDE = []
        _GRIDS_OF_NOISE_AMPLITUDE = []
        
        
        #-----------SCAN SPOT CALCUL---------------------------------------------------
        for scanspot in range(files.shape[0]):  
            
            records = files[scanspot]
            
            print ('--- Computing on scanspot #{} ---'.format(scanspot))
            print (records)
            
            file = '%s/%s.wcp'%(directory,records[0])
            reader = winraw(file)
            reader.parse_header()
    
            nb_sweeps =reader.header['nb_segment'][0]
            sampling = reader.get_signal_sampling_rate()
            
            channel_indexes = np.arange(0,len(reader.header['signal_channels']),1)
            
            raw_sigs = reader.get_analogsignal_chunk(block_index=0,seg_index=1,
                                                     i_start=0,i_stop=-1,
                                                     channel_indexes=channel_indexes)
    
            float_sigs = reader.rescale_signal_raw_to_float(raw_sigs,dtype='float64')
            float_signal = float_sigs[:,ch]
            
            
            points = len(float_signal)
            
            time = np.arange(0,points,1)*1./sampling
            
            if protocol == 'short 32':
                noise_baseline_begin = np.where(time>=0.15)[0][0]
                noise_baseline_end = np.where(time>=0.21)[0][0]
                
                win_baseline_begin = np.where(time>=0.45)[0][0]
                win_baseline_end = np.where(time>=0.5)[0][0]        
                
                window_begin = np.where(time>=0.5)[0][0]
                window_end = np.where(time>=0.7)[0][0]
                
                noise_begin = np.where(time>=0.21)[0][0]
                noise_end = np.where(time>=0.41)[0][0]
    
            if protocol == 'short 96':
                noise_baseline_begin = np.where(time>=0.15)[0][0]
                noise_baseline_end = np.where(time>=0.21)[0][0]
                
                win_baseline_begin = np.where(time>=0.45)[0][0]
                win_baseline_end = np.where(time>=0.5)[0][0]        
                
                window_begin = np.where(time>=0.5)[0][0]
                window_end = np.where(time>=0.7)[0][0]
                
                noise_begin = np.where(time>=0.21)[0][0]
                noise_end = np.where(time>=0.41)[0][0]
                
            elif protocol == 'ultrafast':
                noise_baseline_begin = np.where(time>=0.09)[0][0]
                noise_baseline_end = np.where(time>=0.10)[0][0]
                
                win_baseline_begin = np.where(time>=0.19)[0][0]
                win_baseline_end = np.where(time>=0.20)[0][0]        
                
                window_begin = np.where(time>=0.20)[0][0]
                window_end = np.where(time>=0.30)[0][0]
                
                noise_begin = np.where(time>=0.10)[0][0]
                noise_end = np.where(time>=0.20)[0][0]
                
            elif protocol == 'stim 200':
                noise_baseline_begin = np.where(time>=0.59)[0][0]
                noise_baseline_end = np.where(time>=0.60)[0][0]
                
                win_baseline_begin = np.where(time>=0.19)[0][0]
                win_baseline_end = np.where(time>=0.20)[0][0]        
                
                window_begin = np.where(time>=0.20)[0][0]
                window_end = np.where(time>=0.40)[0][0]
                
                noise_begin = np.where(time>=0.60)[0][0]
                noise_end = np.where(time>=0.80)[0][0]
            
            _SWEEPS = np.zeros((96,len(records),int(window_end-window_begin))) #Remove +1 for paired recordings 
            _NOISE  = np.zeros((96,len(records),int(noise_end-noise_begin)))
        
                   
            _GRID_SIGNAL_AMP = []
            _GRID_NOISE_AMP = []
            
        #---------------MEASURE WINDOW & LEAK REMOVING---------------------------------
            for record in range(len(records)):
                
                file_ = '%s/%s.wcp'%(directory,records[record])
                #Open the file
                r = winraw(file_)
                
                #Parse header = VERY IMPORTANT
                r.parse_header()
            
                for (row,col),site in np.ndenumerate(grid-1):
                    
                    #Get raw binary signals
                    raw_signal = r.get_analogsignal_chunk(block_index=0,seg_index=site,
                                                          i_start=0,i_stop=-1,
                                                          channel_indexes=channel_indexes)
                    
                    #Scale signals to float = apply gain factor
                    signals = r.rescale_signal_raw_to_float(raw_signal,dtype='float64')
                    
                    #Isolate the signal from the cannel of interest
                    signal = signals[:,ch]
                    
                    #Remove leak and append to matrix for later 
                    win_leak = np.mean(signal[win_baseline_begin:win_baseline_end])
                    noise_leak = np.mean(signal[noise_baseline_begin:noise_baseline_end])
                                        
                    noise = signal[noise_begin:noise_end]-noise_leak
        
                    signal_ = signal[window_begin:window_end]-win_leak
               
                    _SWEEPS[site,record,:]=signal_
                    _NOISE[site,record,:]=noise
        
        #-------------------------PLOT AND AMP MEASURES--------------------------------     
        
            if FAST == False : #Do the signals plot
                f2,ax = plt.subplots(gridX,gridY,sharex=True,sharey=True,figsize=(12,12))
                plt.suptitle('%s scanspot %s'%(name,scanspot+1))
                time_ticks = [window_begin,(window_begin+window_end)/2,window_end]
                tick_labels = [0,100,200]
            
            for (row,col),site in np.ndenumerate(grid-1):
                signal_average = np.nanmean(_SWEEPS[site],axis=0)
                noise_average = np.nanmean(_NOISE[site],axis=0)
                
                if FAST == False : #Do the signals plot 
                    ax[row,col].plot(noise_average,color='0.5')        
                    ax[row,col].plot(signal_average,color='blue')
                    ax[row,col].set_title(site+1)
        
                for index in range(len(signal_average)):
                    if signal_average[index] == np.min(signal_average):
                        min_signal = index
                        signal_window = signal_average[min_signal-50:min_signal+50]
        
                        
                for index in range(len(noise_average)):            
                    if noise_average[index] == np.min(noise_average):
                        min_noise = index
                        noise_window = noise_average[min_noise-50:min_noise+50]
                
        
                if len(signal_window) == 100:
                    _SIGNAL_AMP = np.nanmean(signal_window,axis=0)
                else:
                    _SIGNAL_AMP =np.min(signal_average) #BECAUSE MIN AMP IS TO CLOSE FROM ARRAY BORDER AND INDUCES NAN VALUE
                    
                if len(noise_window) == 100:
                    _NOISE_AMP = np.nanmean(noise_window,axis=0)
                else:
                    _NOISE_AMP = np.min(noise_average)
                
                
                if FAST == False : #Do the signals plot            
                    ax[row,col].plot(np.ones(len(noise_average))*_NOISE_AMP,color='black')
                    ax[row,col].plot(np.ones(len(signal_average))*_SIGNAL_AMP,color='red')
        
                _GRID_SIGNAL_AMP.append(_SIGNAL_AMP)
                _GRID_NOISE_AMP.append(_NOISE_AMP)
    
            if FAST == False : #Do the signals plot   
                if SAVING == True:
                    plt.savefig(r'%s/%s_ScanSpot_%s.png'%(savedir,name,scanspot+1))
                else :
                    print ('Fig is not saved')
            
            if CLOSE_FIG == True:
                plt.close()
        #----------------RESHAPE AND CO------------------------------------------------
            _GRID_SIGNAL_AMP = np.asarray(_GRID_SIGNAL_AMP)    
            _GRID_NOISE_AMP = np.asarray(_GRID_NOISE_AMP)
            
            _IMG_GRID_SIGNAL_AMP = np.reshape(_GRID_SIGNAL_AMP,(grid.shape[0],grid.shape[1]))    
            _IMG_GRID_NOISE_AMP = np.reshape(_GRID_NOISE_AMP,(grid.shape[0],grid.shape[1]))
            
            _GRIDS_OF_SIGNAL_AMPLITUDE.append(_IMG_GRID_SIGNAL_AMP)
            _GRIDS_OF_NOISE_AMPLITUDE.append(_IMG_GRID_NOISE_AMP)    
           
        #---------------------GLOBAL RESHAPE & CONCATENATE-----------------------------
        
        if global_reshape.size == 1 : #Only the first scanspot
            
            fake_map = np.empty((grid.shape[0],grid.shape[1]))
            fake_map.fill(np.nan)
            
            _FULL_SIGNAL_MAP = np.concatenate((fake_map,fake_map,_GRIDS_OF_SIGNAL_AMPLITUDE[0],fake_map),axis=1)
        
            _FULL_NOISE_MAP = np.concatenate((fake_map,fake_map,_GRIDS_OF_NOISE_AMPLITUDE[0],fake_map),axis=1)
        
            
          
        elif global_reshape.size == 4: #Full map
        
            _FULL_SIGNAL_MAP = np.concatenate((_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[0]],_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[1]],_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[2]],_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[3]]),axis=1)
        
            _FULL_NOISE_MAP = np.concatenate((_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[0]],_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[1]],_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[2]],_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[3]]),axis=1)
        
        else: #Adapt because of lacking scanspot
            fake_map = np.empty((grid.shape[0],grid.shape[1]))
            fake_map.fill(np.nan)
            
            _FULL_SIGNAL_MAP = np.array([])
            _FULL_NOISE_MAP = np.array([]) 
        
            for i in range(len(global_reshape)):
        
                _ind = global_reshape[i]   
                
                if _ind == global_reshape[0]:        
                    _FULL_SIGNAL_MAP = _GRIDS_OF_SIGNAL_AMPLITUDE[_ind]
                    _FULL_NOISE_MAP = _GRIDS_OF_NOISE_AMPLITUDE[_ind]
        
        
                if _ind != global_reshape[0]:                    
                    _FULL_SIGNAL_MAP = np.concatenate((_FULL_SIGNAL_MAP,_GRIDS_OF_SIGNAL_AMPLITUDE[_ind]),axis=1)
                    _FULL_NOISE_MAP = np.concatenate((_FULL_NOISE_MAP,_GRIDS_OF_NOISE_AMPLITUDE[_ind]),axis=1)
        
                
            if np.array_equal(global_reshape,np.array((1,0))) == True:   #IF SCANSPOT 3 & 4 ARE MISSING IN 2-1-3-4 CONFIG          
                _FULL_SIGNAL_MAP = np.concatenate((_FULL_SIGNAL_MAP,fake_map,fake_map),axis=1)            
                _FULL_NOISE_MAP = np.concatenate((_FULL_NOISE_MAP,fake_map,fake_map),axis=1)  
                
            elif np.array_equal(global_reshape,np.array((1,0,2))) == True:  #IF SCANSPOT 4 is MISSING IN 2-1-3-4 CONFIG          
                _FULL_SIGNAL_MAP = np.concatenate((_FULL_SIGNAL_MAP,fake_map),axis=1)            
                _FULL_NOISE_MAP = np.concatenate((_FULL_NOISE_MAP,fake_map),axis=1)  
                
            elif np.array_equal(global_reshape,np.array((0,1))) == True:  #IF SCANSPOT 3 & 4 are MISSING IN 4-3-1-2 CONFIG          
                _FULL_SIGNAL_MAP = np.concatenate((fake_map,fake_map,_FULL_SIGNAL_MAP),axis=1)            
                _FULL_NOISE_MAP = np.concatenate((fake_map,fake_map,_FULL_NOISE_MAP),axis=1)  
                
            elif np.array_equal(global_reshape,np.array((2,0,1))) == True:  #IF SCANSPOT 4 is MISSING IN 4-3-1-2 CONFIG          
                _FULL_SIGNAL_MAP = np.concatenate((fake_map,_FULL_SIGNAL_MAP),axis=1)            
                _FULL_NOISE_MAP = np.concatenate((fake_map,_FULL_NOISE_MAP),axis=1)  
                
            elif np.array_equal(global_reshape,np.array((0,1,2))) == True:  #IF SCANSPOT 4 is MISSING IN 4-1-2-3 CONFIG          
                _FULL_SIGNAL_MAP = np.concatenate((fake_map,_FULL_SIGNAL_MAP),axis=1)            
                _FULL_NOISE_MAP = np.concatenate((fake_map,_FULL_NOISE_MAP),axis=1)  
        
        
        #---------------------------MAP FLIP ACCORDING TO POSITION-----------------------
                
                
        if orientation== 1:
            print ("Map in position 1: nothing to change")
            
        if orientation== 2:
            #Flip et position----------------------------------------------------------
            _FULL_SIGNAL_MAP = np.fliplr(_FULL_SIGNAL_MAP)	            #Flips map left-Right axis
            _FULL_NOISE_MAP = np.fliplr(_FULL_NOISE_MAP)	            #Flips map left-Right axis
            grid = np.fliplr(grid)				                          #Flips grid on left-right axis
            print ("Map in position 2: Flip on medio-lateral axis")
            
        if orientation==3:
            _FULL_SIGNAL_MAP = np.flipud(_FULL_SIGNAL_MAP)		#Flips map on top-down axis
            _FULL_NOISE_MAP = np.flipud(_FULL_NOISE_MAP)		    #Flips map on top-down axis   
            grid = np.flipud(grid)	
            print ("Map in position 3: Flip on top-down axis")   
                
        if orientation== 4:
            #Flip et position----------------------------------------------------------
            _FULL_SIGNAL_MAP = np.fliplr(_FULL_SIGNAL_MAP)	
            _FULL_SIGNAL_MAP = np.flipud(_FULL_SIGNAL_MAP)
            _FULL_NOISE_MAP = np.fliplr(_FULL_NOISE_MAP)	
            _FULL_NOISE_MAP = np.flipud(_FULL_NOISE_MAP)    
            
            grid = np.fliplr(grid)	
            grid = np.flipud(grid)			#Flips on left-right axis and top bottom axis
            print ("Map in position 4: Flip on both axis")
        
        
        #----------------------NOISE HISTOGRAM & MEAN----------------------------------
        noise = np.ravel(_FULL_NOISE_MAP)
        noise = noise[np.logical_not(np.isnan(noise))]
        
        mean_noise = np.mean(noise)
        
        plt.figure()
        plt.suptitle('%s noise distribution'%name)
        
        reverse_noise = [item*-1 for item in noise if item<=0]
        
        n,bins,patches, = plt.hist(reverse_noise, number_bin_fit,density=1,facecolor='gray',alpha=0.5, label='noise distrubution') #Histogram
        
        #For fit of gaussian noise
        def gaussian(x, mu, sigma):
            return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))
        
        def half_gauss(x,sigma):
            return 2*sigma/np.pi*np.exp(-x**2*sigma**2/np.pi)
                
        x = bins[:-1]
        y = n
        n = len(x)
        
        px = np.linspace(np.min(x),np.max(x),n)
        
        popt, pcov = curve_fit(half_gauss,x,y)
        sigma = popt[0]

        
        py = half_gauss(px,sigma)
        
        plt.plot(px,py, linewidth=2,color='orange',label='half gaussian fit')
        
        spline = UnivariateSpline(px, py-np.max(py)/2,s=0) #Plane for HWHM
        root_x = float(spline.roots()[0])
        root_y = half_gauss(root_x,sigma)
        
        if aPrioriSigmaValue != None:
            root_x = aPrioriSigmaValue
        
        plt.scatter(root_x,root_y,label='HWHM={}'.format(round(root_x,2)),s=30,color='red')
        plt.plot([root_x,root_x],[0,root_y],linestyle='--', color='red')
        plt.plot([0,root_x],[root_y,root_y],linestyle='--', color='red')
        plt.legend(loc='best')
        
        if SAVING == True:
            plt.savefig(r'%s/%s_Charge_Map_Noise_Hist.png'%(savedir,name))
          
        if CLOSE_FIG == True:
            plt.close()
        
        
        #--------------------------------2D ZSCORE-------------------------------------
        sigma_2D = np.full((gridX,frameY),np.abs(root_x))
        noise_2D = np.full((gridX,frameY),np.abs(mean_noise))
        Zscore_2D = (np.abs(_FULL_SIGNAL_MAP)-noise_2D)/sigma_2D
        
        #-----------------------------1D PATTERN & ZSCORE------------------------------
        Amp_pattern = np.max(np.abs(_FULL_SIGNAL_MAP),axis=0)
        
        mean_array = np.ones(Zscore_2D.shape[1])*np.abs(mean_noise)
        sigma_array = np.ones(Zscore_2D.shape[1])*np.abs(root_x)
        Zscore_1D = (Amp_pattern-mean_array)/sigma_array
        
        #Zscore_1D = np.max(Zscore_2D,axis=0)
        #----------------------------------FIGURE PLOT---------------------------------
        
        fig1,ax = plt.subplots(5,1, sharex=False, sharey=False,figsize=(8.27,11.69))		
        					 
        complete_grid = np.concatenate((grid, grid, grid, grid), axis = 1)
        
        fig1.suptitle('%s' %name, fontsize=20)
        
        v = ax[4].imshow(np.abs(_FULL_SIGNAL_MAP),interpolation='nearest', cmap='hot',vmin = np.abs(root_x))
        cbar1 = plt.colorbar(v,orientation='horizontal')
        ax[4].set_title('2D Max Amplitude',loc='right')
        
        mask = Zscore_2D>=ZSOCRE_CUT				#2D Zscore
        ax[3].imshow(mask, interpolation='none',cmap='hot')
        ax[3].set_title('2D Zscore',loc='right')
        
        zcolor = 'green'; zalpha = 0.5
        ax[2].axvspan(zebrin[0],zebrin[1],color=zcolor,alpha=zalpha)
        ax[2].axvspan(zebrin[2],zebrin[3],color=zcolor,alpha=zalpha)
        ax[2].axvspan(zebrin[2],zebrin[3],color=zcolor,alpha=zalpha)
        ax[2].axvspan(zebrin[5],zebrin[6],color=zcolor,alpha=zalpha)
        
        ax[2].axvspan(zebrin[4]-0.2,zebrin[4]+0.2,color='red',alpha=zalpha) #the cell
        
        
        ax[2].set_xlabel('Distance (P1- normalized - cp aligned)')	
        ax[2].set_title('Zebrin Bands',loc='right')
        ax[2].set_xlim(-lb, rb)
        
        
        x = np.arange(0,len(Amp_pattern),1)
        ax[0].fill_between(x,0,Amp_pattern,color='lightblue',label='Amplitude',interpolate=True)
        ax[0].set_title('Amplitude pattern',loc='right')
        ax[0].set_xlim(0,len(Amp_pattern)-1)
        ax[0].set_ylabel('Absolute Max Amplitude (pA)')
        
        ax[1].fill_between(x,0,Zscore_1D,color='0.4',linewidth=2,label='Zscore')
        z_lim = np.ones(len(Zscore_1D))*ZSOCRE_CUT
        ax[1].fill_between(x,z_lim,Zscore_1D,where=Zscore_1D>=z_lim,color='lightcoral',linewidth=2,label='Zscore',interpolate=True)
        ax[1].set_title('1D Zscore',loc='right')
        ax[1].set_xlim(0,len(Zscore_1D)-1)
        ax[1].set_ylabel('Zscore')
        
        if SAVING == True:
            plt.savefig(r'%s/%s_Amp_Map_Zebrin.pdf'%(savedir,name))
        
        if CLOSE_FIG == True:
            plt.close()
        #-----------------------------SAVINGS------------------------------------------
        
        if SAVING == True :
            np.savetxt(r"%s/%s_Amp_zscore_max_OK.csv" %(savedir,name),Zscore_1D, delimiter=',')
            np.savetxt(r"%s/%s_Amp_zscore_2D_OK.csv" %(savedir,name),Zscore_2D, delimiter=',')
            
            np.savetxt(r"%s/%s_Amp_max_OK.csv" %(savedir,name),Amp_pattern, delimiter=',')
            np.savetxt(r"%s/%s_Amp_2D_OK.csv" %(savedir,name),_FULL_SIGNAL_MAP, delimiter=',')
            
            np.savetxt(r"%s/%s_Amp_Noisemap_OK.csv" %(savedir,name),_FULL_NOISE_MAP, delimiter=',')
            np.savetxt(r"%s/%s_Amp_Sigma_OK.csv" %(savedir,name),sigma_2D, delimiter=',')
            
            np.savetxt(r"%s/%s_Zebrin_OK.csv" %(savedir,name),zebrin, delimiter=',')
            
            np.savetxt(r"%s/%s_Positions_OK.csv" %(savedir,name),positions, delimiter=',')
            np.savetxt(r"%s/%s_Positions_cp_centered_OK.csv" %(savedir,name),positions_cp_centered, delimiter=',')
            
            np.savetxt(r"%s/%s_Files_List.csv" %(savedir,name),files,delimiter=',',fmt="%s")
            
        else : 
            print ("No data was saved")
            
        DONE.append(name)
        
    except Exception : 
        import traceback
        traceback.print_exc()
        continue

print ('These maps have been computed successfully')
print (DONE)
print ('{} maps on {} have been computed'.format(len(DONE),len(NAMES)))

print ([x for x in NAMES if x not in DONE], 'These maps were not computed ')