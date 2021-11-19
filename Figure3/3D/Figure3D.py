# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 18:36:02 2020

@author: ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sn
import pingouin as pg
from scipy import stats
from scipy import integrate
import sys 

file = 'C:/Users/klab/Downloads/DataSource_Spaeth_Bahuguna_et_al (2).xlsx'

saveDir = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/REVISION CODE/00_FINAL_CODE/Figure3/3D'
#Redirect output
sys.stdout = open('{}/Local_Statistics.txt'.format(saveDir),'w')

palette = ['royalblue','0.8','lightcoral']

fig, ax = plt.subplots(1,4, sharex=False, sharey=True)

df = pd.read_excel(file, header=1, sheet_name='Figure 3c',na_values='').loc[:23]

#Easy first, peak amplitude at D15
postOp15 = df[['day 15','Condition']]

sn.boxplot(x='Condition',y='day 15',data=df,ax=ax[0],palette=['lightcoral','0.5','royalblue'])
#sn.swarmplot(x='Condition',y='post_op_15',data=postOp15,ax=ax[0],color='black')

#[1] --------------------Stats on PostOp15 peak-------------------------------------
postOp15_KW = stats.kruskal(df.loc[postOp15['Condition']=='cuff']['day 15'].values,
                            df.loc[postOp15['Condition']=='sham']['day 15'].values,
                            df.loc[postOp15['Condition']=='control']['day 15'].values)



print('Analysis of behavioral features')
print ('[1]--------------------- Post Op 15 ------------------------')



print('')

print('Multi condition test')
print ('Kruskal Wallis')
print (postOp15_KW)
if postOp15_KW[1] <= 0.05:
    print ('Post Hoc MWU tests')
    wt = postOp15.loc[postOp15['Condition']=='control']['day 15'].values
    sh = postOp15.loc[postOp15['Condition']=='sham']['day 15'].values
    cu = postOp15.loc[postOp15['Condition']=='cuff']['day 15'].values
    
    wtVSsh = stats.mannwhitneyu(wt, sh)
    print ('CTRL vs SHAM stat={} ;  p-val={}'.format(wtVSsh[0],round(wtVSsh[1],5)))
    
    wtVScu = stats.mannwhitneyu(wt, cu)
    print ('CTRL vs CUFF stat={} ; p-val={}'.format(wtVScu[0],round(wtVScu[1],5)))

    shVScu = stats.mannwhitneyu(sh, cu)
    print ('SHAM vs CUFF stat={} ; p-val={}'.format(shVScu[0],round(shVScu[1],5)))


#[2]---------------------------------AUCs----------------------------------
#Early : bl to day 4 / late : day 4 to 21
print ('')
print ('[2]----------------------- AUCs ----------------------------')

#Iterate on each condition 
AUC_EARLY, AUC_LATE, RECOVERY = [],[],[]
for condition in ['control','sham','cuff']: 
    print ("")
    print (condition, 'mean', 'std')
    early, late, recover = [],[],[]
    subDf = df.loc[df['Condition']==condition]
    print (subDf.mean())
    print (subDf.std())
    
    for row in range(subDf.shape[0]): 
        
        earlyProfile = subDf.iloc[row][['day 2','day 4','day 9']]
        lateProfile = subDf.iloc[row][['day 9','day 14','day 15','day 21']]
        recoveryProfile = subDf.iloc[row][['day 21','day 28','day 33']]
        
        early.append(integrate.trapz(earlyProfile.values))
        late.append(integrate.trapz(lateProfile.values))
        recover.append(integrate.trapz(recoveryProfile.values))


    print ('Avg +/- STD in early, late and recovery phase')
    print (np.mean(early),'+/-',np.std(early))
    print (np.mean(late),'+/-',np.std(late))
    print (np.mean(recover),'+/-',np.std(recover))
        
    AUC_EARLY.append(early)
    AUC_LATE.append(late)
    RECOVERY.append(recover)
    

    
#Plot the distributions 
dfEarly = pd.DataFrame(AUC_EARLY, index=['CTRL','SHAM','CUFF']).T
dfLate = pd.DataFrame(AUC_LATE, ['CTRL','SHAM','CUFF']).T
dfRec = pd.DataFrame(RECOVERY, ['CTRL','SHAM','CUFF']).T

#--------------Save the DFs---------------------------------------------------
dfEarly.to_excel('{}/_AUC_EARLY.xlsx'.format(saveDir))
dfLate.to_excel('{}/_AUC_LATE.xlsx'.format(saveDir))
dfRec.to_excel('{}/_AUC_RECOVERY.xlsx'.format(saveDir))
postOp15.to_excel('{}/_POST_OP15.xlsx'.format(saveDir))


sn.boxplot(data=dfEarly,ax=ax[1],palette=palette)
#sn.swarmplot(data=dfEarly,ax=ax[1],color='black')
ax[1].set_title('AUC Early')

sn.boxplot(data=dfLate,ax=ax[2],palette=palette)
#sn.swarmplot(data=dfLate,ax=ax[2],color='black')
ax[2].set_title('AUC Late')

sn.boxplot(data=dfRec,ax=ax[3],palette=palette)
#sn.swarmplot(data=dfLate,ax=ax[2],color='black')
ax[3].set_title('AUC Recovery')

#Do the stats
print('')
print('NON PARAMETRIC TESTS')
print('Kruskal Wallis AUC ealry')
earlyKW = stats.kruskal(AUC_EARLY[0],AUC_EARLY[1],AUC_EARLY[2])
print (earlyKW)
print ('Post how MWU tests')

for dist1, cond1 in zip(AUC_EARLY,['CTRL','SHAM','CUFF']):
    
    for dist2, cond2 in zip(AUC_EARLY,['CTRL','SHAM','CUFF']):
    
        mwu = stats.mannwhitneyu(dist1, dist2)
        print ('{} vs {} MWU stat={} ;  p-val={}'.format(cond1,cond2,mwu[0],mwu[1]))

print('')
print('Kruskal Wallis AUC late')
lateKW = stats.kruskal(AUC_LATE[0],AUC_LATE[1],AUC_LATE[2])
print(lateKW)
print ('Post hoc mwu')
for dist1, cond1 in zip(AUC_LATE,['CTRL','SHAM','CUFF']):
    
    for dist2, cond2 in zip(AUC_LATE,['CTRL','SHAM','CUFF']):
    
        mwu = stats.mannwhitneyu(dist1, dist2)
        print ('{} vs {} MWU stat={} ;  p-val={}'.format(cond1,cond2,mwu[0],mwu[1]))
        
print('')
print('Kruskal Wallis AUC recovery')
lateKW = stats.kruskal(RECOVERY[0],RECOVERY[1],RECOVERY[2])
print(lateKW)
print ('Post hoc mwu')
for dist1, cond1 in zip(RECOVERY,['CTRL','SHAM','CUFF']):
    
    for dist2, cond2 in zip(RECOVERY,['CTRL','SHAM','CUFF']):
    
        mwu = stats.mannwhitneyu(dist1, dist2)
        print ('{} vs {} MWU stat={} ;  p-val={}'.format(cond1,cond2,mwu[0],mwu[1]))
        
        
print ('')
print ('Levene test on AUC early')
levene_early = stats.levene(AUC_EARLY[0],AUC_EARLY[1],AUC_EARLY[2])
print (levene_early)

print ('')
print ('Levene test on AUC late')
levene_late = stats.levene(AUC_LATE[0],AUC_LATE[1],AUC_LATE[2])
print (levene_late)

print ('')
print ('Levene test on RECOVERY')
levene_rec = stats.levene(RECOVERY[0],RECOVERY[1],RECOVERY[2])
print (levene_rec)

print ('')
print ('Bartlett test on AUC early')
bartlett_early = stats.bartlett(AUC_EARLY[0],AUC_EARLY[1],AUC_EARLY[2])
print (bartlett_early)
        
        
        
        

print('')
print('PARAMETRIC TESTS')
print('Anova AUC ealry')
earlyKW = stats.f_oneway(AUC_EARLY[0],AUC_EARLY[1],AUC_EARLY[2])
print (earlyKW)
print ('Post how t test tests')

for dist1, cond1 in zip(AUC_EARLY,['CTRL','SHAM','CUFF']):
    
    for dist2, cond2 in zip(AUC_EARLY,['CTRL','SHAM','CUFF']):
    
        mwu = stats.ttest_ind(dist1, dist2)
        print ('{} vs {} t test ind , stat= {} ; p-val={}'.format(cond1,cond2,mwu[0],mwu[1]))

print('')
print('Anova AUC late')
lateKW = stats.f_oneway(AUC_LATE[0],AUC_LATE[1],AUC_LATE[2])
print(lateKW)
print ('Post hoc t test')
for dist1, cond1 in zip(AUC_LATE,['CTRL','SHAM','CUFF']):
    
    for dist2, cond2 in zip(AUC_LATE,['CTRL','SHAM','CUFF']):
    
        mwu = stats.ttest_ind(dist1, dist2)
        print ('{} vs {} MWU stat={}, p-val={}'.format(cond1,cond2,mwu[0],mwu[1]))
        

        
print('')
print('Anova AUC Recovery')
lateKW = stats.f_oneway(RECOVERY[0],RECOVERY[1],RECOVERY[2])
print(lateKW)
print ('Post hoc t test')
for dist1, cond1 in zip(RECOVERY,['CTRL','SHAM','CUFF']):
    
    for dist2, cond2 in zip(RECOVERY,['CTRL','SHAM','CUFF']):
    
        mwu = stats.ttest_ind(dist1, dist2)
        print ('{} vs {} MWU stat={}, p-val={}'.format(cond1,cond2,mwu[0],mwu[1]))

        
        

        
        
        
        
        
fig.savefig('{}/PostOp15_AUCearly_late_Recovery.pdf'.format(saveDir))

    




