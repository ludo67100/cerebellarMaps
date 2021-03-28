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
import bct
from collections import Counter 
import matplotlib.cm as cm

import sys
import seaborn as sns
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler,normalize
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.lines as mlines

matplotlib.rcParams['pdf.fonttype'] = 42

sys.path.append("../common/")
import analyze as anal

# Raw data
data_dir = "/home/bahuguna/Work/Isope_data/Isope_data_cerebellar_maps/"
# Processed data 
data_target_dir = "../data/"
fig_target_dir = "../figs/"

electrophys = "ELECTROPHY"
behavior = "BEHAVIOR"
ENR_7 = ["180517","190517","160517","170517","200517","270820","280820"]
ENR_19 = ["220319","230319","210319","171219","090920","100920"]

zone_names = ["B_contra","AX_contra","Alat_contra","Amed_contra","Amed_ipsi","Alat_ipsi","AX_ipsi","B_ipsi"]

sub_ipsi_contra = sys.argv[1]


if sub_ipsi_contra == "n":  # whole graph properties
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_with_behavior_pandas_all.csv")
    graph_props = ["modularity_index","module_degree_zscore","participation_pos","local_assortativity_pos_whole"]

elif sub_ipsi_contra == "y": # graph properties on the resolution of hemispheric level
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_with_behavior_pandas_sub_ipsi_contra_all.csv")
    graph_props = ["modularity_index"] +  [x+"_"+h  for x in ["module_degree_zscore","participation_pos","local_assortativity_pos"] for h in ["ipsi","contra"] ]
elif sub_ipsi_contra == "zone_wise":
    graph_prop_df = pd.read_csv(data_target_dir+"merged_zone_wise_graph_props_avg_seeds.csv")
    graph_props = ["modularity_index"] + [y+"-"+x  for x in ["module_degree_zscore","participation_pos","local_assortativity_pos_whole"] for y in zone_names ]
elif sub_ipsi_contra == "semi_zone_wise": # Instead of values for every zone, only use the mean and variabnce of the zone wise trajectory in order zone_names
    graph_prop_df = pd.read_csv(data_target_dir+"merged_zone_wise_graph_props_avg_seeds.csv")
    graph_props = []
    for x in ["module_degree_zscore","participation_pos","local_assortativity_pos_whole"]:
        temp_prop = [y+"-"+x  for y in zone_names]
        sub_data = graph_prop_df[temp_prop]
        mean = np.mean(np.array(sub_data),axis=1)
        var = np.var(np.array(sub_data),axis=1)
        graph_props.append(x+"-"+"mean")
        graph_prop_df[x+"-"+"mean"] = mean
        graph_props.append(x+"-"+"var")
        graph_prop_df[x+"-"+"var"] = var




sub_data_temp = graph_prop_df[graph_props+["names","subtype","gamma"]]

non_string_columns = graph_props

# Average over animals
sub_data_temp["short_names"] = [ x.split('-')[0]  for x in sub_data_temp["names"]]
temp_ma = dict()
for k in non_string_columns:
    temp_ma[k] = []
sub_data = anal.mean_over_gammas(["names","gamma"],sub_data_temp,temp_ma,non_string_columns)
for k in non_string_columns:
    sub_data[k] = sub_data[k].astype('float')
sub_data["short_names"] = [ x.split('-')[0]  for x in sub_data["names"]]


mean_over_gamma_summary = dict()
temp_gm = dict()
for x in non_string_columns:
    temp_gm[x] = []

sub_data2_gamma = anal.mean_over_gammas("names",sub_data,temp_gm,non_string_columns)

for x in non_string_columns:
    sub_data2_gamma[x] = sub_data2_gamma[x].astype('float')

sub_data2_gamma["short_names"] = [ x.split('-')[0]  for x in sub_data2_gamma["names"]]

# Now zscore
sub_data2_zs_gp = pd.DataFrame()
for x in non_string_columns : 
    sub_data2_zs_gp[x] = (sub_data2_gamma[x]-np.mean(sub_data2_gamma[x]))/np.std(sub_data2_gamma[x])


sub_data2_zs_gp["subtype"] = sub_data2_gamma["subtype"]
sub_data2_zs_gp["short_names"] =  sub_data2_gamma["short_names"]

uniq_subtype = np.unique(sub_data2_zs_gp["subtype"])
target_clus_labels = np.ones(len(sub_data2_zs_gp["subtype"]))

for i,st in enumerate(uniq_subtype):
    ind = np.where(np.array(sub_data2_zs_gp["subtype"])==st)
    target_clus_labels[ind] = i+1

sub_data2_zs_gp["target_clus_labels"] = target_clus_labels
if sub_ipsi_contra == "n":
    tit = "Whole graph properties"
elif sub_ipsi_contra == "y":
    tit = "Ipsi/Contra graph properties"
elif sub_ipsi_contra == "zone_wise":
    tit = "Zone wise graph properties"
elif sub_ipsi_contra == "semi_zone_wise":
    tit = "Zone wise macro properties"


sns.set(font_scale=1.5)
sns.set_style(style='white')
seed = np.random.randint(0,9999999)
sub_data2_zs_gp = sub_data2_zs_gp.dropna()

seed = 9216226
print(seed)
np.random.seed(seed)


for st in ["subtype"]:
        
    np.random.seed(seed)
    tsne = TSNE(n_components=2, verbose=1, perplexity=10, n_iter=3000,init='random',metric='euclidean',learning_rate=150.0,early_exaggeration=24.)
    tsne_results = tsne.fit_transform(sub_data2_zs_gp[non_string_columns])
    df_plot = pd.DataFrame(columns=["tsne-2d-0","tsne-2d-1"]+[st])
    df_plot[st] = np.array(sub_data2_zs_gp[st])
    df_plot["tsne-2d-0"] = tsne_results[:,0]
    df_plot["tsne-2d-1"] = tsne_results[:,1]

    fig = plt.figure(figsize=(16,16))
    fig.suptitle(tit,fontsize=15,fontweight='bold')
    t1 = fig.add_subplot(111)

    uniq_st = np.unique(df_plot[st])
    colors = cm.get_cmap('tab10',9)
    palette = [ colors(i)   for i in np.arange(9)]


    markers=dict()
    mark_list = ['o','s','^','v','P','D','X','*','<','d']
    for i,k in enumerate(uniq_st):
        markers[k] = mark_list[i] 
    g = sns.scatterplot(x="tsne-2d-0", y="tsne-2d-1",hue=st,palette=sns.color_palette(palette[:len(uniq_st)]),data=df_plot,legend="full",ax=t1,s=200,markers=markers,style=st,alpha=0.5,style_order=uniq_st,hue_order=uniq_st)
    g.legend_.remove()
    
    for i,k in enumerate(uniq_st):
        df_plot_st = df_plot.loc[df_plot[st]==k]

        anal.plot_scatter_plot_errorbars(data=df_plot_st,props=["tsne-2d-0","tsne-2d-1"],color=palette[i],marker=markers[k],ax=t1,ms=20)

   
    g.axes.set_ylabel(g.axes.get_ylabel(),fontsize=20,fontweight='bold')
    g.axes.set_xlabel(g.axes.get_xlabel(),fontsize=20,fontweight='bold')
    g.axes.set_xticklabels([])
    g.axes.set_yticklabels([])
   
    legs = [ mlines.Line2D([],[],color=palette[i],marker=markers[k],linestyle='None',markersize=15,label=k) for i,k in enumerate(uniq_st)]

    fig.legend(legs,uniq_st,loc='best',fontsize=20)
    fig.subplots_adjust(left=0.06,bottom=0.06,right=0.96,top=0.93)
    fig.savefig(fig_target_dir+"tsne_all_subtypes_"+sub_ipsi_contra+"_"+st+"_seeds.png")


