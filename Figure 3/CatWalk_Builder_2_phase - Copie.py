# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:38:53 2019
@author: ludov
"""

#------------------------------ZE CODE IS BELOW--------------------------------

#---------------------------------Functions------------------------------------
def remove_leak(signal,bl_start,bl_stop):
    '''
    Removes leak from trace
    signal (1D array) = the data
    bl_start (int) = index start
    bl_stop (int) = index stop
    
    returns signal with bl around 0
    '''
    import numpy as np
    
    leak = np.nanmean(signal[bl_start:bl_stop],axis=0)
    
    return signal-leak


def analyse_cw_data(data_folder='None',date='None',
                    coord_file='None',
                    savedir='None',
                    removeleak=True,
                    end_plot=False,intermediate_plot=False,
                    saveplot=False,savedata=False):
    
    from neo.rawio import WinWcpRawIO as win 
    from matplotlib import pyplot as plt 
    import numpy as np 
    import pandas as pd
    import scipy as sp
    import glob

    #Main plot
    main_fig, ax = plt.subplots(1,1)
    ax.set_title(date)
    
    LOG_VALUES = []
    
    #Iterate each file, to compute ratios
    for file,idx in zip(glob.glob('{}/*.wcp'.format(data_folder)),range(len(glob.glob('{}/*.wcp'.format(data_folder))))):    
        
        #Load WCP file
        reader = win(file)
        reader.parse_header() #HEADER NEEDS TO BE PARSED !!!!!!!!!!!!!!!
        nb_sweeps = reader.header['nb_segment'][0]
        sampling_rate = reader.get_signal_sampling_rate()
        
        #Load coord dataframe
        coordinates = pd.read_excel(coord_file,header=0,index_col=0,sheet_name=date)
    
        animal_name = file[-7:-4]
        if animal_name == 'WT3':
            continue
        
        if 'S' in animal_name:
            animal_name = animal_name[1:]
        
        print ('----------------')
        print ('Animal: {}'.format(animal_name))
    
        AUCs, RATIOs = [],[]
        
        #Iterate along the sweeps 
        for sweep in range(nb_sweeps):   
            
            print('Trial#{}'.format(sweep+1))
            
            #Operant window : extract time borders and remove discarded trials
            try:
                win_on = float(coordinates.loc[animal_name,sweep+1].split(':')[0])/1000
                win_off = float(coordinates.loc[animal_name,sweep+1].split(':')[1])/1000
                
            except ValueError: #Cause 'DISCARDED' string cannot be converted to float
                print ('Trial has been discarded')
                continue
    
            #Get the raw sigs, binary format 
            #i_start, i_stop = index to start, to stop
            raw_sigs = reader.get_analogsignal_chunk(block_index=0, seg_index=sweep, 
                                                     i_start=int(win_on*sampling_rate), 
                                                     i_stop=int(win_off*sampling_rate),
                                                     channel_indexes=[0,1])
            
            #Convert to float64 to rescale and get corret y values for both channels 
            left_sigs = reader.rescale_signal_raw_to_float(raw_sigs, dtype='float64')[:,0].ravel()
            right_sigs = reader.rescale_signal_raw_to_float(raw_sigs, dtype='float64')[:,1].ravel()
            
            #Check both signals for a given sweep have the same size... That's better indeed
            assert len(left_sigs) == len(right_sigs),"left and right signals are not of equal length, that's embarassing"
            
            #Compute the time vector
            time_vector = np.arange(0,len(left_sigs),1)*1./sampling_rate
        
            #Baseline
            bl_start = np.where(time_vector>=0.0)[0][0]
            bl_stop =  np.where(time_vector>=0.004)[0][0]
            
            if removeleak == True:
                #Get the signals, remove leak
                clean_left_signal = remove_leak(left_sigs,bl_start,bl_stop)
                clean_right_signal = remove_leak(right_sigs,bl_start,bl_stop)*-1
                
            else:
                clean_left_signal = left_sigs
                clean_right_signal = right_sigs*-1
                
        
            #Sinusoid =  left_signal+(right_signal)*-1
            sinusoid = clean_left_signal+clean_right_signal
            
            ratio = np.abs(sp.trapz(clean_left_signal,dx=1./sampling_rate))/np.abs(sp.trapz(clean_right_signal,dx=1./sampling_rate))
            RATIOs.append(ratio)
            print ('ratio =',ratio)
            
            #Integer sinusoid
            sin_AUC = sp.trapz(sinusoid,dx=1./sampling_rate)
            AUCs.append(sin_AUC)
            print ('Sinusoid AUC =',sin_AUC)
            
            #The plot
            if intermediate_plot == True:
                plt.figure()
                plt.xlabel('Time(s)')
                plt.ylabel('Pressure (AU)')
                plt.title('{}, {}, Trial#{}, L/R ratio={:.2f}'.format(date,animal_name,sweep+1,ratio))
                plt.plot(time_vector,clean_left_signal,color='0.3',label='Left')
                plt.plot(time_vector,clean_right_signal,color='0.7',label='right')
                plt.plot(time_vector,sinusoid, color='orange',label='sinusoid',linestyle='--',linewidth=1)
                plt.legend(loc='best')
                
                if saveplot == True: 
                    plt.savefig('{}/{}_{}_Trial{}.pdf'.format(savedir,date,animal_name,sweep+1))
                    
        #Log the ratios, average and plot 
        individual_logs = np.log(RATIOs)

        LOG_VALUES.append(individual_logs)
        
        if savedata == True: 
            np.savetxt('{}/{}_{}_LOGS.csv'.format(savedir,date,animal_name),individual_logs)
            np.savetxt('{}/{}_{}_RAW_RATIO.csv'.format(savedir,date,animal_name),RATIOs)
        
        if end_plot == True:
            ax.scatter(np.ones(len(individual_logs))*idx,individual_logs,label=animal_name)
            ax.set_ylabel('Log L/R ratio')
            ax.legend(loc='best')
            
        else: 
            plt.close()
        
    return LOG_VALUES
        
#--------------------------------The Code--------------------------------------
        
if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    import numpy as np
    
    DATA = []
    
    saving = True
    
    fig1, ax = plt.subplots(1)
    
    data_dir = 'U:/RAW DATA/CatWalk/2e phase/Sham_Wt'
    days = ['Baseline1','Baseline2','Baseline3','Post2','Post4','Post9','Post14','Post15','Post21','Post28','Last_Trial']
    dates = ['BL1-13112019','BL2-14112019','BL3-19112019','POST2-23112019','POST4-25112019','POST9-29-11-2019','POST14-06-12-2019','POST15-07-12-2019','POST21-14-12-2019','POST28-21-12-2019','LAST']

#    days = ['Baseline3','Post2','Post4','Post9','Post14','Post15','Post21','Post28','Last_Trial']
#    dates = ['BL3-19112019','POST2-23112019','POST4-25112019','POST9-29-11-2019','POST14-06-12-2019','POST15-07-12-2019','POST21-14-12-2019','POST28-21-12-2019','LAST']


    condition = 'SHAM'

    coord_file = 'U:/01_ANALYSIS/CatWalk/2e phase/windows.xlsx'
    savedir = 'U:/01_ANALYSIS/CatWalk/00_FOR_PAPER'
    
    for day, date in zip(days,dates):
        
        print (day)
    
        data_folder = '{}/{}/{}'.format(data_dir,day,condition)
        date = date

        #If True, will plot individual traces
        intermediate_plot = False
        
        #If True, will save the plots in savedir
        saveplot = False
    
        #If True, will save individual excel sheets with logs 
        savedata = False
        
        ze_data = analyse_cw_data(data_folder=data_folder,
                                  date=date,
                                  coord_file=coord_file,
                                  savedir=savedir,
                                  removeleak=False,
                                  end_plot=True)

        experiment_day = []
        for mouse in ze_data :
            experiment_day.append(np.nanmean(mouse))
            
        DATA.append(experiment_day)
    
    for i in range(np.asarray(DATA).T.shape[0]):
        ax.plot(np.asarray(DATA).T[i],color='0.5')
    
    ax.plot(np.nanmean(DATA,axis=1),label=condition)
    
    if saving == True: 
        
        import pandas as pd 
        
        if condition == 'SHAM':
            
            animal_list = ['S0','S1','S2','S3','S4']

            cuff_df = pd.DataFrame(data=np.asarray(DATA).T,index=np.unique(animal_list),columns=days)
            
            cuff_df.to_excel('{}/SHAM_SECOND_SET_DATA.xlsx'.format(savedir),sheet_name='SHAM_2')
        
    
