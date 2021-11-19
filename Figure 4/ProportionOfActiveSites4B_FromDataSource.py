# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 14:20:00 2021

@author: klab
"""

import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 7})


file = 'D:/000_PAPER/000/FinalResub/Source_Data_Spaeth_Bahuguna_et_al.xlsx'


import pandas as pd 
import seaborn as sn 


df = pd.read_excel(file, header=0, sheet_name='Figure 4b')


sn.boxplot(data=df, x='Condition',y='Proportion of Active sites (%)', showfliers=False, order=['control','short training','long training','early sham','early cuff', 'adapted sham','adapted cuff'])
sn.swarmplot(data=df, x='Condition',y='Proportion of Active sites (%)',size=8,color='0.2',order=['control','short training','long training','early sham','early cuff', 'adapted sham','adapted cuff'])

import pingouin as pg 


kruskal = pg.kruskal(data=df, between='Condition', dv = 'Proportion of Active sites (%)')
print (kruskal)