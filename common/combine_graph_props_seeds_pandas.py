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


#  
data_target_dir = "data/"
fig_target_dir = "figs/"

data_type = sys.argv[1]

print(data_type)

if data_type == "subtype":

    electrophys = "ELECTROPHY"
    # Raw data
    data_dir = "For Paper/EPHYS/Adaptive_Dataset/"

    subtypes = os.listdir(data_dir)
    #data_2d = pickle.load(open(data_target_dir+"data_2d_maps.pickle","rb"))
    #data = pd.read_csv(data_target_dir+"meta_data.csv")
    
    files = glob.glob(data_target_dir+"graph_properties_norm_*.pickle")

elif data_type == "development":
    development = "DEVELOPMENT"
    # Raw data
    data_dir = "For Paper/EPHYS/Development_Dataset/"
    subtypes = os.listdir(data_dir) # Just the name of the variable is subtypes, its actually days
    #data_2d = pickle.load(open(data_target_dir+"data_2d_maps_days.pickle","rb"))
    #data = pd.read_csv(data_target_dir+"meta_data_days.csv")


    files = glob.glob(data_target_dir+"graph_properties_days_norm_*.pickle")


num_or_size = "num" # num of clusters or size of the largest cluster
gamma_re_arrange = 0.34



gammas = np.arange(0.0,1.5,0.17)
cmaps = [cm.get_cmap('Reds',len(gammas)+10), cm.get_cmap('Blues',len(gammas)+10), cm.get_cmap('Greens',len(gammas)+10), cm.get_cmap('Purples',len(gammas)+10),cm.get_cmap('Greys',len(gammas)+4),cm.get_cmap('pink_r',len(gammas)+10)]


graph_prop_simps = dict()

percentile = 70

dat_type = data_type
print(files)
print(len(files))
for f in files:
    seed = f.split('/')[-1].split('_')[-1].split('.')[0]
    graph_properties = pickle.load(open(f,"rb"))

    graph_prop_df = pd.DataFrame(columns=["modularity_index","gamma","participation_pos","participation_neg","local_assortativity_pos_whole","module_degree_zscore","total_amplitude","average_amplitude","percentage_active_sites","names"]+[dat_type])


    temp_dict = dict()
    for x in list(graph_prop_df.keys()):
        temp_dict[x] = []

    for i,st in enumerate(subtypes):
        st_list_cov=[]
        st_mods_list_cov=[]
        st_list_corr=[]
        st_mods_list_corr=[]
        norms =[]
        tot_amp=[]
        avg_amp = []
        per_act_sit = []
        graph_prop_simps[st] = dict()
        participation_pos = []
        participation_neg = []
        loc_ass_pos = []
        #loc_ass_neg = []
        zscore = []
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
                elif "total_amplitude" in y1:
                    tot_amp.append(graph_properties[st]["modularity"][x]["total_amplitude"])
                elif "average_amplitude" in y1:
                    avg_amp.append(graph_properties[st]["modularity"][x]["average_amplitude"])
                elif "percentage_active_sites" in y1:
                    per_act_sit.append(graph_properties[st]["modularity"][x]["percentage_active_sites"])
                elif "participation" in y1 and "whole" in y1:
                    participation_pos.append(graph_properties[st]["modularity"][x]["participation_whole"][0])
                    participation_neg.append(graph_properties[st]["modularity"][x]["participation_whole"][1])
                elif "zscore" in y1 and "whole" in y1:
                    zscore.append(graph_properties[st]["modularity"][x]["module_degree_zscore_whole"])
                elif "local" in y1:
                    loc_ass_pos.append(graph_properties[st]["modularity"][x]["local_assortativity_whole"])
                elif y1 == "cov" or y1 == "corr":
                    mod_indices = graph_properties[st]["modularity"][x][y1][0]
                    num_mods = [len(y) for y in  graph_properties[st]["modularity"][x][y1][1]]
                    # If num_mods are zero just go to next data point, because if this empty, causes problems, while slicing by gammas
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
                    

        graph_prop_simps[st]["participation_pos"] = participation_pos
        graph_prop_simps[st]["participation_neg"] = participation_neg
        graph_prop_simps[st]["module_degree_zscore"] = zscore

        print(len(norms),len(st_list_corr))

        nz_inds = np.unique(nz_inds)
        if len(norms) > len(st_list_corr):
            graph_prop_simps[st]["st_list_corr_norm"] = np.array(norms)[nz_inds]
            graph_prop_simps[st]["total_amplitude"] = np.array(tot_amp)[nz_inds]
            graph_prop_simps[st]["average_amplitude"] = np.array(avg_amp)[nz_inds]
            graph_prop_simps[st]["percentage_active_sites"] = np.array(per_act_sit)[nz_inds]

        else:
            graph_prop_simps[st]["st_list_corr_norm"] = np.array(norms)
            graph_prop_simps[st]["total_amplitude"] = np.array(tot_amp)
            graph_prop_simps[st]["average_amplitude"] = np.array(avg_amp)
            graph_prop_simps[st]["percentage_active_sites"] = np.array(per_act_sit)

        if len(loc_ass_pos) > len(st_list_corr):
            graph_prop_simps[st]["local_assortativity_pos_whole"] = np.array(loc_ass_pos)[nz_inds]
        else:
            graph_prop_simps[st]["local_assortativity_pos_whole"] = np.array(loc_ass_pos)


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
                temp_dict["total_amplitude"].append(np.array(tot_amp)[nz_inds])
                temp_dict["average_amplitude"].append(np.array(avg_amp)[nz_inds])
                temp_dict["percentage_active_sites"].append(np.array(per_act_sit)[nz_inds])
                temp_dict["participation_pos"].append(np.array(graph_prop_simps[st]["participation_pos"])[nz_inds,k])
                temp_dict["participation_neg"].append(np.array(graph_prop_simps[st]["participation_neg"])[nz_inds,k])
                temp_dict["module_degree_zscore"].append(np.array(graph_prop_simps[st]["module_degree_zscore"])[nz_inds,k])

            else:
                temp_dict["total_amplitude"].append(np.array(tot_amp))
                temp_dict["average_amplitude"].append(np.array(avg_amp))
                temp_dict["percentage_active_sites"].append(np.array(per_act_sit))
                temp_dict["participation_pos"].append(np.array(graph_prop_simps[st]["participation_pos"])[:,k])
                temp_dict["participation_neg"].append(np.array(graph_prop_simps[st]["participation_neg"])[:,k])
                temp_dict["module_degree_zscore"].append(np.array(graph_prop_simps[st]["module_degree_zscore"])[:,k])

            
            if len(names) > len(st_list_corr):
                temp_dict["names"].append(np.array(graph_prop_simps[st]["names"])[nz_inds])
            else:
                temp_dict["names"].append(np.array(graph_prop_simps[st]["names"]))
            
            temp_dict["local_assortativity_pos_whole"].append(np.array(graph_prop_simps[st]["local_assortativity_pos_whole"]))

            count+=len(np.array(st_list_corr)[:,:,k][:,0]) 
        
        temp_dict[dat_type].append( [st for i3 in np.arange(0,count)])

        print(st)
        print(len(st_list_cov))
        print(len(st_list_corr))


    graph_prop_df["modularity_index"] = np.hstack(temp_dict["modularity_index"])
    graph_prop_df["total_amplitude"] = np.hstack(temp_dict["total_amplitude"])
    graph_prop_df["average_amplitude"] = np.hstack(temp_dict["average_amplitude"])
    graph_prop_df["percentage_active_sites"] = np.hstack(temp_dict["percentage_active_sites"])

    graph_prop_df["participation_pos"] = np.hstack(temp_dict["participation_pos"])
    graph_prop_df["local_assortativity_pos_whole"] = np.hstack(temp_dict["local_assortativity_pos_whole"])
    graph_prop_df["participation_neg"] = np.hstack(temp_dict["participation_neg"])
    graph_prop_df["module_degree_zscore"] = np.hstack(temp_dict["module_degree_zscore"])
    graph_prop_df["names"] = np.hstack(temp_dict["names"])


    graph_prop_df["gamma"] = np.hstack(temp_dict["gamma"])

    graph_prop_df[dat_type] = np.hstack(temp_dict[dat_type])

    graph_prop_df = graph_prop_df.replace([np.inf, -np.inf], np.nan)

    if data_type == "subtype":
        graph_prop_df.to_csv(data_target_dir+"graph_properties_pandas_for_behav_"+seed+".csv")
    elif data_type == "development":
        graph_prop_df.to_csv(data_target_dir+"graph_properties_pandas_for_behav_days_"+seed+".csv")

    graph_prop_df_nonan = graph_prop_df.dropna(axis=0)

    if data_type == "subtype":
        graph_prop_df_nonan.to_csv(data_target_dir+"graph_properties_pandas_"+seed+".csv")
    elif data_type == "development":
        graph_prop_df_nonan.to_csv(data_target_dir+"graph_properties_pandas_days_"+seed+".csv")



if data_type == "subtype":
    files1 = glob.glob(data_target_dir+"graph_properties_pandas_for_behav_[0-9]*.csv")
    files2 = glob.glob(data_target_dir+"graph_properties_pandas_[0-9]*.csv")
elif data_type == "development":
    files1 = glob.glob(data_target_dir+"graph_properties_pandas_for_behav_days_[0-9]*.csv")
    files2 = glob.glob(data_target_dir+"graph_properties_pandas_days_[0-9]*.csv")

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
    merge_df1.to_csv(data_target_dir+"graph_properties_pandas_for_behav_all.csv") # everything

    merge_df2 = merge_df_seeds(files2)
    merge_df2.to_csv(data_target_dir+"graph_properties_pandas_all.csv") # nonan

elif data_type == "development":
    merge_df1 = merge_df_seeds(files1)
    merge_df1.to_csv(data_target_dir+"graph_properties_pandas_for_behav_days_all.csv")

    merge_df2 = merge_df_seeds(files2)
    merge_df2.to_csv(data_target_dir+"graph_properties_pandas_days_all.csv")




if data_type == "subtype":
    post_fix = ""
elif data_type == "development":
    post_fix = "_days_"




