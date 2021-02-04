# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:01:06 2019

Gets significant amplitudes per maps, and average it according to A, AX, X and B zones defined by Jan Voogd
UPDATE : INCLUDES INCOMPLETE MAPS, replace missing values with NANs

Amplitudes are centered to the average


@author: Ludovic.SPAETH
"""


import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


#-----------------------DIRECTORY & INFO---------------------------------------
#--------------CHECK ALL THESE INFO BEFORE RUNNING----------------------------- 
# FILES AND DIRECTORY ------------------------------------------------------------------------------------------------------------

SAVING = True

groups = ['WT','ENR1','ENR2','LC','LS','EC','ES']
colors = ['skyblue','limegreen','green','lightcoral','black','orange','purple']
sheets = ['WT','ENR','ENR','LC','LS','EC','ES']


inputDir = 'E:/000_PAPER/Amplitude_Analysis/Sigma'

savedir = 'E:/000_PAPER/Amplitude_Analysis/Sigma/04_ZONES_CORRELATION/Norm Max'

z_score_limit = 3.09


def AppendMean(OriginList,TargetList):

    #Function that appends mean of list if elements in it or 0 if empty
    if len(OriginList)==0:
        TargetList.append(0)
    else :  
        if np.nan in OriginList : 
            TargetList.append(np.nan)
        else:
            TargetList.append(np.nanmean(OriginList))


def micron_to_norm(measure, P1_size):   
    return float(measure)*100./P1_size


for group,sheet in zip(range(len(groups)),sheets):

    condition = groups[group]
    df  = pd.read_excel('E:/000_PAPER/Mesures_ZII_HighRes_WT_ENR_EC_LC_ES_LS.xlsx',
                        sheet_name=sheet,index_col=0,header=1)  #DataFrame containing zebrin band file 
    
    directory = '{}/{}'.format(inputDir,condition)
    
    #files = glob.glob(r'%s/*_Amp_2D_OK.csv'%directory)  #Directory
    names = os.listdir(directory)
               
    List_bands = ['Contra_1','Contra_2', 'Contra_3', 'Contra_4','Ipsi_1', 'Ipsi_2', 'Ipsi_3', 'Ipsi_4']
    
    _MEAN =[]
    _MED = []
    _MAX = []
    _COUNT = []
    _SUM = []
    _LAT = []
    _TOT_AMP= []
    
    _names_OK= [] #to append names of maps that are able to run into analysis 
    
    #All amps 
    A_med_i, A_lat_i, Ax_i, B_i, A_med_c, A_lat_c, Ax_c, B_c = [],[],[],[],[],[],[],[]  

    #Significant sites only        
    sA_med_i, sA_lat_i, sAx_i, sB_i, sA_med_c, sA_lat_c, sAx_c, sB_c = [],[],[],[],[],[],[],[]  
    
    sC4,sC3,sC2,sC1,sI1,sI2,sI3,sI4 = [],[],[],[],[],[],[],[]  
    
    cA_med_i, cA_lat_i, cAx_i, cB_i, cA_med_c, cA_lat_c, cAx_c, cB_c = [],[],[],[],[],[],[],[] #For count active vs inactive
     
    for name in range(len(names)):  
        
        subDir = '{}/{}'.format(directory,names[name])
        
        Amp2D = np.abs(np.genfromtxt('{}/{}_Amp_2D_OK.csv'.format(subDir,names[name]),delimiter=','))
        Zscore2D = np.abs(np.genfromtxt('{}/{}_Amp_zscore_2D_OK.csv'.format(subDir,names[name]),delimiter=','))
        z_score = np.abs(np.genfromtxt('{}/{}_Amp_zscore_max_OK.csv'.format(subDir,names[name]),delimiter=','))
        positions = np.genfromtxt('{}/{}_Positions_OK.csv'.format(subDir,names[name]),delimiter=',')        
        
#        Amp2D = np.genfromtxt(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\AMP_ANALYSIS\001_FOR_PAPER\INDIVIDUAL_MAPS\%s\%s_Amp_2D_OK.csv"%(condition,names[name]),delimiter=',')
#        Zscore2D = np.genfromtxt(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\AMP_ANALYSIS\001_FOR_PAPER\INDIVIDUAL_MAPS\%s\%s_Amp_zscore_2D_OK.csv"%(condition,names[name]),delimiter=',')
#        z_score = np.genfromtxt(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\AMP_ANALYSIS\001_FOR_PAPER\INDIVIDUAL_MAPS\%s\%s_Amp_max_OK.csv"%(condition,names[name]),delimiter=',')
#        positions = np.genfromtxt(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\AMP_ANALYSIS\001_FOR_PAPER\INDIVIDUAL_MAPS\%s\%s_Positions_OK.csv"%(condition,names[name]),delimiter=',')

##### NO NORMALIZATION = Comment everything         
        
##### Normed max
        MAX_AMP = np.nanmax(Amp2D)
        Amp2D = Amp2D / MAX_AMP

###### CENTERED:
#        AVERAGE_AMP = np.nanmean(Amp2D)
#        Amp2D = Amp2D - AVERAGE_AMP
#        
###### REDUIT:        
#        STD_AMP = np.nanstd(Amp2D)
#        Amp2D = Amp2D / STD_AMP
#        
        
        #Previous way
        #Amp2D = Amp2D / float(np.min(Amp2D.ravel()))
        
#        #Norme centre
#        normed_Amp2D = Amp2D/np.nanmax(Amp2D)
#        
#        
#        norm_avg = np.nanmean(normed_Amp2D)
#        print (norm_avg)
#        
#        del Amp2D
#        
#        Amp2D = normed_Amp2D-norm_avg
        
#################################################

        #Getting the bands
        _BANDS_norm =  df.loc[['%s norm_P1-'%names[name]],'P2- contra':'P2- ipsi']  #Selects values in xcl file from P2- contra to P2- ipsi
        P1mAVERAGE = 290.5 #Average size of P1-
        
        
        
        B_contra = -233.
        Ax_contra = -133.
        Alat_contra = -108.
        Amed_contra = -58.
        
        Amed_ipsi = 50.
        Alat_ipsi = 100.
        Ax_ipsi = 125.
        B_ipsi = 225
            


        _names_OK.append(names[name])
        
        
        incontra4,incontra3,incontra2,incontra1,inipsi1,inipsi2,inipsi3,inipsi4 = [],[],[],[],[],[],[],[]
        
        print ('-----{}-----'.format(names[name]))        
        for position, idx in zip(positions, range(len(positions))):
            #P2- contra
            if B_contra <= position < Ax_contra :
                for i in range(len(Amp2D[:,idx])):                    
                    if np.isnan(Amp2D[i,idx]) == True:
                        incontra4.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                incontra4.append(Amp2D[i,idx])       
                            
            #P2+ contra
            if Ax_contra <= position < Alat_contra :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        incontra3.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                incontra3.append(Amp2D[i,idx])  
                            
            #P1- contra
            if Alat_contra <= position < Amed_contra :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        incontra2.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                incontra2.append(Amp2D[i,idx])  
            
            #P1p
            if Amed_contra <= position < 0 :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        incontra1.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                incontra1.append(Amp2D[i,idx])   
            

            #Local Cluster
            if 0 <= position < Amed_ipsi :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        inipsi1.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                inipsi1.append(Amp2D[i,idx])  
                            
            #P1- ipsi
            if Amed_ipsi <= position < Alat_ipsi :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        inipsi2.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                inipsi2.append(Amp2D[i,idx])  

            #P2+ ipsi
            if Alat_ipsi <= position < Ax_ipsi :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        inipsi3.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                inipsi3.append(Amp2D[i,idx])  
                            
            #P2- ipsi
            if Ax_ipsi <= position < B_ipsi :
                for i in range(len(Amp2D[:,idx])):
                    if np.isnan(Amp2D[i,idx]) == True:
                        inipsi4.append(np.nan)
                    
                    else:
                            if Zscore2D[i,idx] >= z_score_limit:
                                inipsi4.append(Amp2D[i,idx])  
                            
                            


                            
                            
        #Append mean of significant sites in Lists
        AppendMean(incontra4,sC4)
        AppendMean(incontra3,sC3)
        AppendMean(incontra2,sC2)
        AppendMean(incontra1,sC1)            
        AppendMean(inipsi1,sI1)
        AppendMean(inipsi2,sI2)
        AppendMean(inipsi3,sI3)
        AppendMean(inipsi4,sI4)
        
    
    Data = np.vstack((sC4,sC3,sC2,sC1,sI1,sI2,sI3,sI4)).transpose()
    cols = ['B_contra','Ax_Contra','A_lat_contra','A_med_contra',
            'A_med_ipsi','A_lat_ipsi','Ax_ipsi','B_ipsi']

    #--------------------------------SAVE-------------------------------------------------    
    
    # Create dataframe to hold data
    Output = pd.DataFrame(Data,index=_names_OK,columns=cols)
        
    if SAVING == True :

        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter('%s/%s_Map_2D_Average_Amp_per_zones_wNANS_%s.xlsx'%(savedir,condition,z_score_limit), engine='xlsxwriter',options={'nan_inf_to_errors': True})
        
        #_Output.style.applymap(hightlight_cols,subset=pd.IndexSlice[:,['P1mC_count','P1mI_count','P2mC_count','P1mC_mean','P1mI_mean','Loc_Clust_Mean','P1mC_med','P1mI_med','P2mI_med','P1mI_sum','P1mC_sum']])
        Output.to_excel(writer, sheet_name='%s_Average_Amplitudes'%condition,na_rep = 'nan') 
        writer.save()
        
    else :
        print ('No datasheet has been saved')
                            
        
        

                    
    

