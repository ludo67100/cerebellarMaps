import os
import glob
import numpy as np
import pylab as pl
import scipy.io as sio
# for_Jyotika.m
from copy import copy, deepcopy
import pickle
import matplotlib.cm as cm
import pdb
import h5py
import pandas as pd
import scipy.stats as sp_st
import sys
import seaborn as sns


# Raw data
data_dir = "/home/bahuguna/Work/Isope_data/Isope_data_cerebellar_maps/"

data_target_dir = "../data/"
fig_target_dir = "../figs/"

electrophys = "ELECTROPHY"
behavior = "BEHAVIOR"

sub_ipsi_contra = sys.argv[1]

behavior_enrichment = pd.read_excel(data_dir+behavior+"/"+"Enrichment.xlsx")
gammas = np.round(np.arange(0.0,1.5,0.17),2)

day_label_order = list(behavior_enrichment.keys())[1:]

enrichment_df = pd.DataFrame(columns=["mouse","time",'Distance','intercept','slope','maximum_distance','total_distance','average_distance',"total_days","short-names"])

temp_df = dict()
for k in list(enrichment_df.keys()):
    temp_df[k] = []

days = behavior_enrichment.keys()[1:]
for i in np.arange(len(behavior_enrichment)):
    x = behavior_enrichment.iloc[i]
    for d in days:
        temp_df["Distance"].append(float(x[d]))
        temp_df["time"].append(d)

    y_dist = np.array(np.array(x)[1:]).astype('float')
    ind_nonan = np.where(np.isnan(y_dist)==False)[0]
    y_dist1 = y_dist[ind_nonan]
    x_days = np.arange(0,len(y_dist1))

    coef = np.polyfit(x_days,y_dist1,1)
    max_dist = np.max(y_dist1)
    tot_dist = np.sum(y_dist1)


    temp_df["mouse"].append([np.array(x["Mouse"])  for i in np.arange(len(days)) ])
    temp_df["short-names"].append([np.array(x["Mouse"].split('_')[1])  for i in np.arange(len(days)) ])
    temp_df["intercept"].append([ coef[1] for i in np.arange(len(days)) ])
    temp_df["slope"].append([ coef[0] for i in np.arange(len(days)) ])
    temp_df["maximum_distance"].append([ max_dist for i in np.arange(len(days)) ])
    temp_df["total_distance"].append([ tot_dist for i in np.arange(len(days)) ])
    temp_df["average_distance"].append([ tot_dist/len(y_dist1) for i in np.arange(len(days)) ])
    temp_df["total_days"].append([len(y_dist1) for i in np.arange(len(days)) ])


for k in list(enrichment_df):
    enrichment_df[k] = np.hstack(temp_df[k])



enrichment_df.to_csv(data_target_dir+"Enrichment_df.csv")

fig = pl.figure(figsize=(16,16))
t1 = fig.add_subplot(111)
g1 = sns.lineplot(x='time',y='Distance',hue='mouse',data=enrichment_df,linewidth=2.5,palette='nipy_spectral',marker='o',ax=t1,sort=False)
fig.savefig(fig_target_dir+"Enrichment_distances.png")

if sub_ipsi_contra == "n":
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_pandas_for_behav_all.csv")
else:
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_all.csv")

graph_prop_enr = graph_prop_df.loc[graph_prop_df["subtype"]=="ENR"]
graph_prop_enr["short-names"] = [ x.split('-')[0] for x in np.array(graph_prop_enr["names"]) ]

graph_prop_enr_behav = pd.merge(enrichment_df,graph_prop_enr,right_on='short-names',left_on='short-names')


if sub_ipsi_contra == "n":
    graph_prop_enr_behav.to_csv(data_target_dir+"graph_properties_behavior_enr_all.csv")
else:
    graph_prop_enr_behav.to_csv(data_target_dir+"graph_properties_behavior_enr_sub_ipsi_contra_all.csv")


