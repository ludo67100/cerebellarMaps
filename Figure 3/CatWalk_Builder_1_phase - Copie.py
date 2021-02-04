# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 13:03:27 2019

@author: Ludovic.SPAETH
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 12:53:22 2019

@author: Ludovic.SPAETH
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 16:20:20 2018

@author: ludovic.spaeth
"""

import numpy as np
from matplotlib import pyplot as plt 
from scipy.integrate import trapz
import matplotlib 
matplotlib.rcParams['pdf.fonttype'] = 42
from scipy import stats

animal_list = ['AI','AII','AIII','BI','BII','BIII','BIV',
               'CI','CII','CIII','CIV','DI','DII','DIII',
               'DIV'] #List your animals here
               
conditions =['BL_Day_8','Cuff_Day_2','Cuff_Day_4','Cuff_Day_9','Cuff_Day_14','Cuff_Day_15','Cuff_Day_21','Cuff_Day_28','Cuff_Day_33']


is_sham = [True, False, False, False, False, False, False,
           True, True, True, True, False, False, False,
           False, True] #True if animal is sham, False if animal if cuffed


directory = 'U:/01_ANALYSIS/CatWalk/Manip/1 mois/Individual Analysis'

savedir = 'U:/01_ANALYSIS/CatWalk/00_FOR_PAPER'

saving = True 

CUFF = np.zeros(len(conditions))
CUFF_SEM = np.zeros(len(conditions))
CUFF_SD = np.zeros(len(conditions))
SHAM = np.zeros(len(conditions))
SHAM_SEM = np.zeros(len(conditions))
SHAM_SD = np.zeros(len(conditions))

sham_idx, cuff_idx = [],[]


for cond in range(len(conditions)):
    
    cuff_values = []
    sham_values = []

#    plt.figure()
#    plt.title(conditions[cond])
#    plt.plot(np.zeros(len(animal_list)),color='black',linestyle='--')
#    plt.ylim(-5,5)    
    
    data = np.genfromtxt('%s/%s/%s_LIST_RATIOS.csv'%(directory,conditions[cond],conditions[cond]),delimiter=',')
    
    log_data = np.log(data)
    
#    for animal in range(len(animal_list)):
#
#        if is_sham[animal] == True:
#            color='green'
#        else:
#            color='orange'
#        
#        plt.scatter(animal,np.nansum((log_data[animal])),color=color)
        
    for animal in range(len(animal_list)):
        if is_sham[animal]==True:
            sham_values.append(np.nanmean((log_data[animal])))
            sham_idx.append(animal_list[animal])
        else:
            cuff_values.append(np.nanmean((log_data[animal])))
            cuff_idx.append(animal_list[animal])
        
    cuff_values = np.asarray(cuff_values)    
    sham_values = np.asarray(sham_values)       
    


    test = stats.mannwhitneyu(cuff_values,sham_values)
    print (conditions[cond],test[1])
    

    cuff_values=cuff_values.flatten()
    sham_values=sham_values.flatten()
    
    if cond == 0: #Its the first loop : 
        CUFF_TABLE = cuff_values
        SHAM_TABLE = sham_values

    else: 
        CUFF_TABLE = np.vstack((CUFF_TABLE,cuff_values))
        SHAM_TABLE = np.vstack((SHAM_TABLE,sham_values))
        
    CUFF[cond] = np.nanmean(cuff_values)
    SHAM[cond] = np.nanmean(sham_values)
    
    
    
    
    CUFF_SEM[cond] = np.nanstd(cuff_values)/np.sqrt(len(cuff_values))
    SHAM_SEM[cond] = np.nanstd(sham_values)/np.sqrt(len(sham_values))
    
    
    CUFF_SD[cond] = np.nanstd(cuff_values)
    SHAM_SD[cond] = np.nanstd(sham_values)
    
    
#Transpose cuff and sham table
CUFF_TABLE = CUFF_TABLE.transpose()
SHAM_TABLE = SHAM_TABLE.transpose()


plt.figure()

for i in range(CUFF_TABLE.shape[0]): 
    plt.plot(CUFF_TABLE[i],color = 'orange')
    
for i in range(SHAM_TABLE.shape[0]):
    plt.plot(SHAM_TABLE[i],color ='green')


plt.plot(np.zeros(len(conditions)),color='black',linestyle='--')

x = np.arange(len(conditions))
plt.plot(CUFF,color='orange',linewidth=1)
plt.fill_between(x,CUFF-CUFF_SEM,CUFF+CUFF_SEM,alpha=0.2,color='orange')

plt.plot(SHAM,color='green',linewidth=1)
plt.fill_between(x,SHAM-SHAM_SEM,SHAM+SHAM_SEM,alpha=0.2,color='green')

plt.xticks(x,conditions,rotation='45')
plt.ylabel('AVG log left_right ratio +/- SEM')
plt.xlabel('time (days)')
#
#np.savetxt('U:/01_ANALYSIS/CatWalk/Manip/1 mois/00_AVG_GROUP/CUFF_AVG.csv',CUFF,delimiter=',')
#np.savetxt('U:/01_ANALYSIS/CatWalk/Manip/1 mois/00_AVG_GROUP/SHAM_AVG.csv',SHAM,delimiter=',')
#
#np.savetxt('U:/01_ANALYSIS/CatWalk/Manip/1 mois/00_AVG_GROUP/CUFF_SD.csv',CUFF_SD,delimiter=',')
#np.savetxt('U:/01_ANALYSIS/CatWalk/Manip/1 mois/00_AVG_GROUP/SHAM_SD.csv',SHAM_SD,delimiter=',')

if saving == True : 
    import pandas as pd 
    
    #For cuff mice
    cuff_df = pd.DataFrame(data=CUFF_TABLE,index=np.unique(cuff_idx),columns=conditions)
    
    cuff_df.to_excel('{}/CUFF_FIRST_SET_DATA.xlsx'.format(savedir),sheet_name='CUFF_1')

    #For sham mice
    sham_df = pd.DataFrame(data=SHAM_TABLE,index=np.unique(sham_idx),columns=conditions)
    
    sham_df.to_excel('{}/SHAM_FIRST_SET_DATA.xlsx'.format(savedir),sheet_name='SHAM_1')
    
    








