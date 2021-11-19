# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 18:58:46 2021

@author: klab
"""

import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 7})
import pandas as pd 
import seaborn as sn 
import matplotlib.pyplot as plt 
import pingouin as pg 
import sys 
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

file = 'C:/Users/klab/Downloads/DataSource_Spaeth_Bahuguna_et_al (2).xlsx'
saveDir = 'D:/000_PAPER/00_ANSWER_TO_REVIEWERS/REVISION CODE/00_FINAL_CODE/Figure2/2HIJK'

fig, ax  = plt.subplots(1,4,figsize=(10,2.5))

sys.stdout = open('{}/Local_Statistics.txt'.format(saveDir),'w')

order = ['P9P10','P12P13','P14P18','P30P40']
colors= ['mediumorchid','lightsteelblue','dodgerblue','navy']
variables = ['modularity_index','participation_pos','module_degree_zscore','local_assortativity_pos_whole']

graphDf = pd.read_excel(file, header=1,sheet_name='Figure 2h-i-j-k')


#Open lists to store the results of the stat comparison
kruskalComp, kruskalLabels = [],[]
mwuComp, mwuLabelsA,mwuLabelsB, mwuVariable = [],[],[],[]
descriptions = []

for i in range(len(variables)):
    
    print()
    print('----------'+variables[i]+'----------')

    subDf = graphDf[[variables[i],'development']]
    #subDf[variables[i]] = 100*subDf[variables[i]]
    
    
    sn.boxplot(y=variables[i],x='development',data=subDf,order=order,ax=ax[i],palette=colors)
    sn.swarmplot(y=variables[i],x='development',data=subDf,order=order,color='black',ax=ax[i])
    
    
    kruskalWallis = pg.kruskal(data=subDf, dv=variables[i], between='development')
    print(kruskalWallis)
    
    kruskalLabels.append('{}'.format(variables[i]))
    kruskalComp.append(kruskalWallis)
    
    #Post hoc 2 by 2 comparisons with MWU
    for distA in order:
        
        ref = subDf.loc[subDf['development']==distA][variables[i]].values
        
        #Get description of the distributions 
        dfStat = subDf.loc[subDf['development']==distA].describe()
        dfStat['Condition'] = [distA for x in range(dfStat.shape[0])]
        descriptions.append(dfStat)
        
        
        for distB in order:
            
            mwuLabelsA.append(distA)
            mwuLabelsB.append(distB)
            mwuVariable.append(variables[i])
            
            compare = subDf.loc[subDf['development']==distB][variables[i]].values
            
            mwu = pg.mwu(ref, compare,tail='two-sided')
            
            mwuComp.append(mwu)
            
            print('{} vs {}'.format(distA, distB))
            print(mwu)
            print()
            
    
    print()
    print()
    
    
#OutPutDf 
kruskalDf = pd.concat(kruskalComp)
kruskalDf['Variable'] = kruskalLabels

mwuDf = pd.concat(mwuComp)
mwuDf['Dist. A'] = mwuLabelsA; mwuDf['Dist. B'] = mwuLabelsB
mwuDf['Variable'] = mwuVariable
mwuDf=mwuDf.replace('P9P10','PND9-10')
mwuDf=mwuDf.replace('P12P13','PND12-13')
mwuDf=mwuDf.replace('P14P18','PND14-18')
mwuDf=mwuDf.replace('P30P40','PND>30')

descriptionDf = pd.concat(descriptions)


#Export dFs

kruskalDf.to_excel('{}/GraphPropDevKruskalWallisTests.xlsx'.format(saveDir))
mwuDf.to_excel('{}/GraphPropDevMannWhitneyUTests.xlsx'.format(saveDir))
descriptionDf.to_excel('{}/GraphPropDevDescription.xlsx'.format(saveDir))
                


plt.tight_layout()
plt.savefig('{}/GraphPropDevBoxplots.pdf'.format(saveDir))