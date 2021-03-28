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
import analyze as anal
import sys
import seaborn as sns
import scipy.stats as sp_stats

# Raw data here
data_dir = "/home/bahuguna/Work/Isope_data/Isope_data_cerebellar_maps/"

data_target_dir = "../data/"
fig_target_dir = "../figs/"

gammas = np.arange(0.0,1.5,0.17)

gamma = sys.argv[1] # "all"(pool over all gammas) or "mean" (average over all gammas) or 0.0 (a particular gamma value)
which_feat = sys.argv[2] # "modularity_index", "participation_pos", "module_degree_zscore", "local_assortativity_pos_whole"
data_type = sys.argv[3]  # "development" or "subtype" (Development data or adaptation data) 
ENR_7 = ["180517","190517","160517","170517","200517","270820","280820"]
ENR_19 = ["220319","230319","210319","171219","090920","100920"]

panel_name = dict({"modularity_index":"H","participation_pos":"I","module_degree_zscore":"J","local_assortativity_pos_whole":"K"})


if data_type == "subtype":
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_pandas_all.csv")
    short_names = [ x.split('-')[0]  for x in graph_prop_df["names"]]
    graph_prop_df["short_names"] = short_names
    phenotypes = np.array(graph_prop_df["subtype"])
    ind_7 = [i    for i,(x,y) in enumerate(zip(phenotypes,np.array(graph_prop_df["short_names"]))) if x == "ENR" and y in ENR_7 ]
    phenotypes[ind_7] = "ENR_7"
    ind_19 = [i    for i,(x,y) in enumerate(zip(phenotypes,np.array(graph_prop_df["short_names"]))) if x == "ENR" and y in ENR_19 ] 
    phenotypes[ind_19] = "ENR_19"

    graph_prop_df["subtype"] = phenotypes
    fig_target_dir = "../figs/"
elif data_type == "development":
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_pandas_days_all.csv")
    fig_target_dir = "../figs/"

if data_type == "subtype":
    colors = ["skyblue","forestgreen","salmon","gray","orange","purple"]
    st_order = ["WT","ENR_7","ENR_19","LC","LS","EC","ES"]
elif data_type == "development":
    colors = ["skyblue","cadetblue","dodgerblue","mediumblue"]
    st_order = ["P9P10","P12P13","P14P18","P30P40"]


#--------------------------------------------------------------------------#
# Average graph properties over gammas and cells of the same animal
#--------------------------------------------------------------------------#

def avg_over_gammas_and_animals(field,data,non_string_columns,metric='mean'):                        
    temp_dat1_wgp = dict()
    for k in non_string_columns:                                                                     
        temp_dat1_wgp[k] = []
        
    if metric == "mean":
        sub_dat1_wgp = anal.mean_over_gammas(field,data,temp_dat1_wgp,non_string_columns,data_type)            
    elif metric == "median":
        sub_dat1_wgp = anal.median_over_gammas(field,data,temp_dat1_wgp,non_string_columns,data_type)          
    for k in non_string_columns:
        sub_dat1_wgp[k] = sub_dat1_wgp[k].astype('float')                                            
        
    return sub_dat1_wgp                                                          





#--------------------------------------------------------------------------#
# Plot the distributions of the graph properties  
#--------------------------------------------------------------------------#


def plot_graph_prop(data,gamma,which_feat,data_type,with_phenotypes=False):
    if type(gamma) == str and gamma == 'all':
        # Pool the results over all gammas
        sub_data = data
        tit = "Pooled over values of gamma"
        figname = "Graph_props_pooled_over_gammas_"+which_feat+"_"+data_type
    elif type(gamma) == str and gamma == "mean":
        tit = "Averaged over values of gamma"
        #figname = "Graph_props_averaged_over_gammas_"+which_feat+"_"+data_type
        figname = "Figure2_Panel"+panel_name[which_feat]+"_"+data_type
        unnamed = [x  for x in data.keys() if "Unnamed" in x]
        non_string_columns = list(set(data.keys())-set(unnamed+ ['names']+[data_type]))
        sub_data = avg_over_gammas_and_animals(["names"],data,non_string_columns)

    else: #type(gamma) == float:
        sub_data = data.loc[data["gamma"]==float(gamma)]
        tit = r'$\gamma$'+" = "+str(gamma)
        figname = "Graph_props_gamma_"+str(gamma)+"_"+which_feat+"_"+data_type

    fig = pl.figure(figsize=(16,16))
    t1 = fig.add_subplot(111)
    g1 = sns.swarmplot(x=data_type,y=which_feat,data=sub_data,ax=t1,order=st_order,color="k",alpha=0.3,size=13)
    g1 = sns.boxplot(x=data_type,y=which_feat,data=sub_data,ax=t1,linewidth=4.5,order=st_order,palette=colors)

    print(which_feat)
    print(gamma)
    for grp in sub_data.groupby(data_type):
        print("===========================================================")
        print(grp[0])
        print("Mean stats")
        print(str(np.mean(grp[1][which_feat]))+"+"+str(np.std(grp[1][which_feat])))
        print("-------------------------------------------------------") 
        print("Median stats")
        print(str(np.median(grp[1][which_feat]))+"+"+str(sp_stats.median_abs_deviation(grp[1][which_feat])))

    
    data_to_save = sub_data[["names","modularity_index","participation_pos","module_degree_zscore","local_assortativity_pos_whole"]+[data_type]]
    data_to_save.to_csv(data_target_dir+"Graph_props_"+which_feat+"_"+data_type+".csv") 

    num_dat_type = len(np.unique(sub_data[data_type]))
    uniq_dat_type = np.unique(sub_data[data_type])

    all_data = []
    greater_sig = dict()
    for i,x1 in enumerate(uniq_dat_type):
        greater_sig[x1] = dict()
        dat1 = np.array(sub_data.loc[sub_data[data_type]==x1][which_feat])
        all_data.append(dat1)
        for j,x2 in enumerate(uniq_dat_type):
            greater_sig[x1][x2] = (0,0)


    feats = []
    
    positions = g1.axes.get_xticks()
    pos_labels = [ x.get_text() for x in g1.axes.get_xticklabels()] 
    for i,x1 in enumerate(uniq_dat_type):
        for j,x2 in enumerate(uniq_dat_type):
            if i == j:
                continue
            feats.append((x1,x2))
            dat1 = np.array(sub_data.loc[sub_data[data_type]==x1][which_feat])
            dat2 = np.array(sub_data.loc[sub_data[data_type]==x2][which_feat])
            t_stat,p_val = sp_stats.mannwhitneyu(dat1,dat2) # Does not assume equal variance
            greater_sig[x1][x2] = (t_stat,p_val)
            if p_val/2. < 0.05 :
                if p_val/2. < 0.01:
                    if p_val < 0.0001:
                        displaystring = r'***'
                    elif p_val < 0.001:
                        displaystring = r'**'
                    else:
                        displaystring = r'*'

                    max_dat = np.max(np.hstack((dat1,dat2)))
                    y_sig = max_dat+(i+j)*0.1*max_dat
                    ind1 = np.where(np.array(pos_labels)==x1)[0][0]
                    ind2 = np.where(np.array(pos_labels)==x2)[0][0]
                    if data_type == "subtype":
                        if x1 == "WT" or x2 == "WT":
                            anal.significance_bar(positions[ind1],positions[ind2],y_sig,displaystring,t1)
                    else:

                        anal.significance_bar(positions[ind1],positions[ind2],y_sig,displaystring,t1)
    # Save it first
    filename = "t_test_"+data_type+"_"+which_feat+"_"+str(gamma)
    pickle.dump(greater_sig,open(data_target_dir+filename+".pickle","wb"))



    g1.axes.set_title(tit,fontsize=15,fontweight='bold')
    g1.axes.set_ylabel(g1.axes.get_ylabel(),fontsize=15,fontweight='bold')
    for x in g1.axes.get_xticklabels():
        x.set_fontsize(20)
        x.set_fontweight('bold')

    g1.axes.set_xlabel("")
    g1.axes.set_ylabel(g1.axes.get_ylabel(),fontsize=20,fontweight='bold')
    #t1.hlines(y=med,xmin=t1.get_xlim()[0],xmax=t1.get_xlim()[1],linestyles='dashed',color='r',linewidth=5.0)

    g1.figure.savefig(fig_target_dir+figname+".png")

    return greater_sig





t_test = plot_graph_prop(graph_prop_df,gamma,which_feat,data_type,False)
print(t_test)

    
    


