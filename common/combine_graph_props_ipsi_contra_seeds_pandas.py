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


# Raw data
data_dir = "/home/bahuguna/Work/Isope_data/Isope_data_cerebellar_maps/"
#  
data_target_dir = "../data/"

fig_target_dir = "../figs/"

num_or_size = "size" # num of clusters or size of the largest cluster
data_type = sys.argv[1]


gamma_re_arrange = 0.34

if data_type == "subtype":

    electrophys = "ELECTROPHY"
    subtypes = os.listdir(data_dir+electrophys)
    data_2d = pickle.load(open(data_target_dir+"data_2d_maps.pickle","rb"))
    data = pd.read_csv(data_target_dir+"meta_data.csv")

    files = glob.glob(data_target_dir+"graph_properties_norm_*.pickle")

elif data_type == "development":

    #fig_target_dir = "/home/bahuguna/Work/Isope_data/Isope_data_analysis/figs/sub_ipsi_contra_development/"
    development = "DEVELOPMENT"
    subtypes = os.listdir(data_dir+development) # Just the name of the variable is subtypes, its actually days
    data_2d = pickle.load(open(data_target_dir+"data_2d_maps_days.pickle","rb"))
    data = pd.read_csv(data_target_dir+"meta_data_days.csv")


    files = glob.glob(data_target_dir+"graph_properties_days_norm_*.pickle")


gammas = np.arange(0.0,1.5,0.17)

graph_prop_simps = dict()

percentile = 70
dat_type = data_type

for f in files:
    seed = f.split('/')[-1].split('_')[-1].split('.')[0]
    graph_properties = pickle.load(open(f,"rb"))

    graph_prop_df = pd.DataFrame(columns=["modularity_index","gamma","norms","participation_pos_whole","participation_pos_ipsi","participation_pos_contra","participation_neg_whole","participation_neg_ipsi","participation_neg_contra","local_assortativity_pos_whole","local_assortativity_pos_ipsi","local_assortativity_pos_contra","module_degree_zscore_whole","module_degree_zscore_ipsi","module_degree_zscore_contra","names"]+[dat_type])

    temp_dict = dict()
    for x in list(graph_prop_df.keys()):
        temp_dict[x] = []
    for i,st in enumerate(subtypes):
        st_list_cov=[]
        st_mods_list_cov=[]
        st_list_corr=[]
        st_mods_list_corr=[]
        norms =[]
        graph_prop_simps[st] = dict()
        
        participation_pos_whole = []
        participation_pos_ipsi = []
        participation_pos_contra = []
        participation_neg_whole = []
        participation_neg_ipsi = []
        participation_neg_contra = []
        loc_ass_pos_whole = []
        loc_ass_pos_ipsi = []
        loc_ass_pos_contra = []
        zscore_whole = []
        zscore_ipsi = []
        zscore_contra = []
        names=[]
        nz_inds = []
        count = 0
        print("==================================================================")
        print(st)
        print("==================================================================")
        for j,x in enumerate(list(graph_properties[st]["modularity"].keys())):
            ind = graph_properties[st]["indices"]
            for y1 in list(graph_properties[st]["modularity"][x].keys()):
                if "norm" in y1:        
                    norms.append(graph_properties[st]["modularity"][x]["norm"])
                elif "participation" in y1 and "whole" in y1:
                    participation_pos_whole.append(graph_properties[st]["modularity"][x]["participation_whole"][0])
                    participation_neg_whole.append(graph_properties[st]["modularity"][x]["participation_whole"][1])
                elif "participation" in y1 and "ipsi" in y1:
                    participation_pos_ipsi.append(graph_properties[st]["modularity"][x]["participation_ipsi"][0])
                    participation_neg_ipsi.append(graph_properties[st]["modularity"][x]["participation_ipsi"][1])
                elif "participation" in y1 and "contra" in y1:
                    participation_pos_contra.append(graph_properties[st]["modularity"][x]["participation_contra"][0])
                    participation_neg_contra.append(graph_properties[st]["modularity"][x]["participation_contra"][1])
                elif "zscore" in y1 and "whole" in y1:
                    zscore_whole.append(graph_properties[st]["modularity"][x]["module_degree_zscore_whole"])
                elif "zscore" in y1 and "ipsi" in y1:
                    zscore_ipsi.append(graph_properties[st]["modularity"][x]["module_degree_zscore_ipsi"])
                elif "zscore" in y1 and "contra" in y1:
                    zscore_contra.append(graph_properties[st]["modularity"][x]["module_degree_zscore_contra"])
              
                elif "local" in y1:
                    loc_ass_pos_whole.append(graph_properties[st]["modularity"][x]["local_assortativity_whole"])
                    loc_ass_pos_ipsi.append(graph_properties[st]["modularity"][x]["local_assortativity_ipsi"])
                    loc_ass_pos_contra.append(graph_properties[st]["modularity"][x]["local_assortativity_contra"])
                elif y1 == "cov" or y1 == "corr":
                    mod_indices = graph_properties[st]["modularity"][x][y1][0]
                    num_mods = [len(y) for y in  graph_properties[st]["modularity"][x][y1][1]]
                    if num_mods[0] == 0:
                        continue
                    num_mods_size = [np.max(y) for y in  graph_properties[st]["modularity"][x][y1][1] if len(y) > 0] 
                    
                    num_mods_greater_size = [ len(np.where(np.array(y) >= np.percentile(y,percentile))[0])  for y in  graph_properties[st]["modularity"][x][y1][1] if len(y) > 0]
                    nz_inds.append(x)
            
                    print(mod_indices)
                    print(num_mods)
                    if "cov" in y1:
                        st_list_cov.append((mod_indices,num_mods,num_mods_size,num_mods_greater_size))
                        st_mods_list_cov.append(graph_properties[st]["modularity"][x][y1][1])
                    elif "corr" in y1:
                        st_list_corr.append((mod_indices,num_mods,num_mods_size,num_mods_greater_size))
                        st_mods_list_corr.append(graph_properties[st]["modularity"][x][y1][1])


        graph_prop_simps[st]["participation_pos_whole"] = participation_pos_whole
        graph_prop_simps[st]["participation_pos_ipsi"] = participation_pos_ipsi
        graph_prop_simps[st]["participation_pos_contra"] = participation_pos_contra
        graph_prop_simps[st]["participation_neg_whole"] = participation_neg_whole
        graph_prop_simps[st]["participation_neg_ipsi"] = participation_neg_ipsi
        graph_prop_simps[st]["participation_neg_contra"] = participation_neg_contra
        graph_prop_simps[st]["module_degree_zscore_whole"] = zscore_whole
        graph_prop_simps[st]["module_degree_zscore_ipsi"] = zscore_ipsi
        graph_prop_simps[st]["module_degree_zscore_contra"] = zscore_contra

        print(len(norms),len(st_list_corr))

        nz_inds = np.unique(nz_inds)
        if len(loc_ass_pos_whole) > len(st_list_corr):
            graph_prop_simps[st]["local_assortativity_pos_whole"] = np.array(loc_ass_pos_whole)[nz_inds]
            graph_prop_simps[st]["local_assortativity_pos_ipsi"] = np.array(loc_ass_pos_ipsi)[nz_inds]
            graph_prop_simps[st]["local_assortativity_pos_contra"] = np.array(loc_ass_pos_contra)[nz_inds]
        else:
            graph_prop_simps[st]["local_assortativity_pos_whole"] = np.array(loc_ass_pos_whole)
            graph_prop_simps[st]["local_assortativity_pos_ipsi"] = np.array(loc_ass_pos_ipsi)
            graph_prop_simps[st]["local_assortativity_pos_contra"] = np.array(loc_ass_pos_contra)

        if len(graph_properties[st]['names']) > len(st_list_corr):
            graph_prop_simps[st]["names"] = np.array(graph_properties[st]['names'])[nz_inds]
        else:
            graph_prop_simps[st]["names"] = np.array(graph_properties[st]['names'])

        if num_or_size == "num":
            ind_prop = 1
        elif num_or_size == "size":
            ind_prop = 2
        for k in np.arange(0,len(gammas)):
            
            temp_dict["modularity_index"].append(np.array(st_list_corr)[:,:,k][:,0])

            nz_inds = np.unique(nz_inds)
            temp_dict["gamma"].append([ np.round(gammas[k],2) for i2 in np.arange(0,len(np.array(st_list_corr)[:,:,k][:,0]))])
            if len(norms) > len(st_list_corr):
                temp_dict["norms"].append(np.array(norms)[nz_inds])
                temp_dict["participation_pos_whole"].append(np.array(graph_prop_simps[st]["participation_pos_whole"])[nz_inds,k])
                temp_dict["participation_pos_ipsi"].append(np.array(graph_prop_simps[st]["participation_pos_ipsi"])[nz_inds,k])
                temp_dict["participation_pos_contra"].append(np.array(graph_prop_simps[st]["participation_pos_contra"])[nz_inds,k])
                temp_dict["participation_neg_whole"].append(np.array(graph_prop_simps[st]["participation_neg_whole"])[nz_inds,k])
                temp_dict["participation_neg_ipsi"].append(np.array(graph_prop_simps[st]["participation_neg_ipsi"])[nz_inds,k])
                temp_dict["participation_neg_contra"].append(np.array(graph_prop_simps[st]["participation_neg_contra"])[nz_inds,k])
                temp_dict["module_degree_zscore_whole"].append(np.array(graph_prop_simps[st]["module_degree_zscore_whole"])[nz_inds,k])
                temp_dict["module_degree_zscore_ipsi"].append(np.array(graph_prop_simps[st]["module_degree_zscore_ipsi"])[nz_inds,k])
                temp_dict["module_degree_zscore_contra"].append(np.array(graph_prop_simps[st]["module_degree_zscore_contra"])[nz_inds,k])

            else:
                temp_dict["norms"].append(np.array(norms))
                temp_dict["participation_pos_whole"].append(np.array(graph_prop_simps[st]["participation_pos_whole"])[:,k])
                temp_dict["participation_pos_ipsi"].append(np.array(graph_prop_simps[st]["participation_pos_ipsi"])[:,k])
                temp_dict["participation_pos_contra"].append(np.array(graph_prop_simps[st]["participation_pos_contra"])[:,k])
                temp_dict["participation_neg_whole"].append(np.array(graph_prop_simps[st]["participation_neg_whole"])[:,k])
                temp_dict["participation_neg_ipsi"].append(np.array(graph_prop_simps[st]["participation_neg_ipsi"])[:,k])
                temp_dict["participation_neg_contra"].append(np.array(graph_prop_simps[st]["participation_neg_contra"])[:,k])
                temp_dict["module_degree_zscore_whole"].append(np.array(graph_prop_simps[st]["module_degree_zscore_whole"])[:,k])
                temp_dict["module_degree_zscore_ipsi"].append(np.array(graph_prop_simps[st]["module_degree_zscore_ipsi"])[:,k])
                temp_dict["module_degree_zscore_contra"].append(np.array(graph_prop_simps[st]["module_degree_zscore_contra"])[:,k])

            if len(names) > len(st_list_corr):
                temp_dict["names"].append(np.array(graph_prop_simps[st]["names"])[nz_inds])
            else:
                temp_dict["names"].append(np.array(graph_prop_simps[st]["names"]))
            
            temp_dict["local_assortativity_pos_whole"].append(np.array(graph_prop_simps[st]["local_assortativity_pos_whole"]))
            temp_dict["local_assortativity_pos_ipsi"].append(np.array(graph_prop_simps[st]["local_assortativity_pos_ipsi"]))
            temp_dict["local_assortativity_pos_contra"].append(np.array(graph_prop_simps[st]["local_assortativity_pos_contra"]))

            count+=len(np.array(st_list_corr)[:,:,k][:,0]) 
            
        temp_dict[dat_type].append( [st for i3 in np.arange(0,count)])

        print(st)
        print(len(st_list_cov))
        print(len(st_list_corr))




    graph_prop_df["modularity_index"] = np.hstack(temp_dict["modularity_index"])
    graph_prop_df["norms"] = np.hstack(temp_dict["norms"])
    graph_prop_df["participation_pos_whole"] = np.hstack(temp_dict["participation_pos_whole"])
    graph_prop_df["participation_pos_ipsi"] = np.hstack(temp_dict["participation_pos_ipsi"])
    graph_prop_df["participation_pos_contra"] = np.hstack(temp_dict["participation_pos_contra"])
    graph_prop_df["local_assortativity_pos_whole"] = np.hstack(temp_dict["local_assortativity_pos_whole"])
    graph_prop_df["local_assortativity_pos_ipsi"] = np.hstack(temp_dict["local_assortativity_pos_ipsi"])
    graph_prop_df["local_assortativity_pos_contra"] = np.hstack(temp_dict["local_assortativity_pos_contra"])
    graph_prop_df["participation_neg_whole"] = np.hstack(temp_dict["participation_neg_whole"])
    graph_prop_df["participation_neg_ipsi"] = np.hstack(temp_dict["participation_neg_ipsi"])
    graph_prop_df["participation_neg_contra"] = np.hstack(temp_dict["participation_neg_contra"])
    graph_prop_df["module_degree_zscore_whole"] = np.hstack(temp_dict["module_degree_zscore_whole"])
    graph_prop_df["module_degree_zscore_ipsi"] = np.hstack(temp_dict["module_degree_zscore_ipsi"])
    graph_prop_df["module_degree_zscore_contra"] = np.hstack(temp_dict["module_degree_zscore_contra"])
    graph_prop_df["names"] = np.hstack(temp_dict["names"])


    graph_prop_df["gamma"] = np.hstack(temp_dict["gamma"])
    graph_prop_df[dat_type] = np.hstack(temp_dict[dat_type])


    graph_prop_df = graph_prop_df.replace([np.inf, -np.inf], np.nan)

    if data_type == "subtype":
        graph_prop_df.to_csv(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_"+seed+".csv")
    elif data_type == "development":
        graph_prop_df.to_csv(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_"+seed+".csv")
    graph_prop_df_nonan = graph_prop_df.dropna(axis=0)


    if data_type == "subtype":
        graph_prop_df_nonan.to_csv(data_target_dir+"graph_properties_pandas_sub_contra_ipsi_"+seed+".csv")
    elif data_type == "development":
        graph_prop_df_nonan.to_csv(data_target_dir+"graph_properties_pandas_sub_contra_ipsi_days_"+seed+".csv")

if data_type == "subtype":
    files1 = glob.glob(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_[0-9]*.csv")
    files2 = glob.glob(data_target_dir+"graph_properties_pandas_sub_contra_ipsi_[0-9]*.csv")
elif data_type == "development":
    files1 = glob.glob(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_[0-9]*.csv")
    files2 = glob.glob(data_target_dir+"graph_properties_pandas_sub_contra_ipsi_days_[0-9]*.csv")

def merge_df_seeds(files):

    for i,f in enumerate(files):
        temp_df = pd.read_csv(f)
        seed = f.split('/')[-1].split('_')[-1].split('.')[0]
        temp_df["seed"] = seed

        if i == 0:
            merge_df = temp_df
        else:
            merge_df = merge_df.append(temp_df)
    return merge_df

if data_type == "subtype":
    merge_df1 = merge_df_seeds(files1)
    merge_df1.to_csv(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_all.csv")

    merge_df2 = merge_df_seeds(files2)
    merge_df2.to_csv(data_target_dir+"graph_properties_pandas_sub_contra_ipsi_all.csv")

elif data_type == "development":
    merge_df1 = merge_df_seeds(files1)
    merge_df1.to_csv(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_all.csv")

    merge_df2 = merge_df_seeds(files2)
    merge_df2.to_csv(data_target_dir+"graph_properties_pandas_sub_contra_ipsi_days_all.csv")



