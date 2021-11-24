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
data_target_dir = "../data/"
fig_target_dir = "../figs/"

data_type = sys.argv[1]

print(data_type)

if data_type == "subtype":

    #fig_target_dir = "/home/bahuguna/Work/Isope_data/Isope_data_analysis/figs/whole_graph_properties/"
    electrophys = "ELECTROPHY"
    # Raw data
    data_dir = "../../COMPLETE_DATASET/For\\ Paper/EPHYS/Adaptive_Dataset/"

    subtypes = os.listdir(data_dir)
    data_2d = pickle.load(open(data_target_dir+"data_2d_maps.pickle","rb"))
    data = pd.read_csv(data_target_dir+"meta_data.csv")
    
    files = glob.glob(data_target_dir+"graph_properties_norm_*.pickle")
    cov_2d_dict =  pickle.load(open(data_target_dir+"covariance_maps_norm.pickle","rb"))

elif data_type == "development":
    development = "DEVELOPMENT"
    # Raw data
    data_dir = "../../COMPLETE_DATASET/For\\ Paper/EPHYS/Development_Dataset/"
    subtypes = os.listdir(data_dir) # Just the name of the variable is subtypes, its actually days
    data_2d = pickle.load(open(data_target_dir+"data_2d_maps_days.pickle","rb"))
    data = pd.read_csv(data_target_dir+"meta_data_days.csv")

    cov_2d_dict =  pickle.load(open(data_target_dir+"covariance_maps_days_norm.pickle","rb"))
    #graph_properties = pickle.load(open(data_dir+"graph_properties_days_norm.pickle","rb"))

    files = glob.glob(data_target_dir+"graph_properties_days_norm_*.pickle")


num_or_size = "num" # num of clusters or size of the largest cluster
gamma_re_arrange = 0.34



gammas = np.arange(0.0,1.5,0.17)
cmaps = [cm.get_cmap('Reds',len(gammas)+10), cm.get_cmap('Blues',len(gammas)+10), cm.get_cmap('Greens',len(gammas)+10), cm.get_cmap('Purples',len(gammas)+10),cm.get_cmap('Greys',len(gammas)+4),cm.get_cmap('pink_r',len(gammas)+10)]


graph_prop_simps = dict()

percentile = 70

dat_type = data_type

for f in files:
    seed = f.split('/')[-1].split('_')[-1].split('.')[0]
    graph_properties = pickle.load(open(f,"rb"))

    graph_prop_df = pd.DataFrame(columns=["modularity_index","gamma","participation_pos","participation_neg","local_assortativity_pos_whole","module_degree_zscore","total_amplitude","average_amplitude","percentage_active_sites","names"]+[dat_type])


    temp_dict = dict()
    for x in list(graph_prop_df.keys()):
        temp_dict[x] = []

    '''
    if data_type == "subtype":
        dat_type = subtypes
    elif data_type == "development":
        dat_type = days
    '''

    for i,st in enumerate(subtypes):
        st_list_cov=[]
        st_mods_list_cov=[]
        st_list_corr=[]
        st_mods_list_corr=[]
        #norms =[]
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
        '''
        for g in gammas:
            within_distances[g] = []
            between_distances[g] = []
        ''' 
        nz_inds = []
        count = 0
        print("==================================================================")
        print(st)
        print("==================================================================")
        for j,x in enumerate(list(graph_properties[st]["modularity"].keys())):
            ind = graph_properties[st]["indices"]
            #names_rn_cn = [  rn+"-"+cn   for rn,cn in zip(np.array(graph_properties[st]["names"][j])[:,0],np.array(graph_properties[st]["names"][j])[:,1]) ]
            #temp_name = np.unique(np.array(graph_properties[st]["names"][j])[:,0])[0]
            #names.append(np.unique(names_rn_cn))
            for y1 in list(graph_properties[st]["modularity"][x].keys()):
                #if "norm" in y1:        
                #    norms.append(graph_properties[st]["modularity"][x]["norm"])
                if "total_amplitude" in y1:
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
                    #loc_ass_neg.append(graph_properties[st]["modularity"][x]["local_assortativity"][1])
                #elif "rich" in y1:
                #    rich_club.append(graph_properties[st]["modularity"][x]["rich_club"])
                elif y1 == "cov" or y1 == "corr":
                    #if data_type == "development":
                    #    pdb.set_trace()
                    mod_indices = graph_properties[st]["modularity"][x][y1][0]
                    #if num_or_size == "num":
                    num_mods = [len(y) for y in  graph_properties[st]["modularity"][x][y1][1]]
                    # If num_mods are zero just go to next data point, because if this empty, causes problems, while slicing by gammas
                    if num_mods[0] == 0:
                        continue
                    #elif num_or_size == "size":
                    '''
                    num_mods_size = [np.max(y) for y in  graph_properties[st]["modularity"][x][y1][1] if len(y) > 0] 
                    
                    num_mods_greater_size = [ len(np.where(np.array(y) >= np.percentile(y,percentile))[0])  for y in  graph_properties[st]["modularity"][x][y1][1] if len(y) > 0]
                    '''
                    nz_inds.append(x)
            
                    print(mod_indices)
                    print(num_mods)
                    #num_mods = len(np.array(graph_properties[st]["modularity"][x][1])[:,0])
                    if "cov" in y1:
                        st_list_cov.append((mod_indices,num_mods,num_mods_size,num_mods_greater_size))
                        st_mods_list_cov.append(graph_properties[st]["modularity"][x][y1][1])
                        #t1.plot(mod_indices,num_mods,'.-',color=cmap1(i+j),alpha=0.5,markersize=8)
                    elif "corr" in y1:
                        st_list_corr.append((mod_indices,num_mods,num_mods_size,num_mods_greater_size))
                        st_mods_list_corr.append(graph_properties[st]["modularity"][x][y1][1])
            #if st == 'ES':
            #    pdb.set_trace()

        graph_prop_simps[st]["st_list_cov"] = st_list_cov
        graph_prop_simps[st]["st_list_corr"] = st_list_corr
        graph_prop_simps[st]["st_mods_list_corr"] = st_mods_list_corr
        graph_prop_simps[st]["st_mods_list_cov"] = st_mods_list_cov
        # To maintain the dimension of (num_subjects, len_gammas)
        #graph_prop_simps[st]["within_module_dists"] =[ [ within_distances[g][ind] for g in gammas ] for ind in np.arange(len(graph_prop_simps[st]["st_mods_list_corr"])) ]
        #graph_prop_simps[st]["between_module_dists"] =[ [ between_distances[g][ind] for g in gammas ] for ind in np.arange(len(graph_prop_simps[st]["st_mods_list_corr"])) ]

        graph_prop_simps[st]["participation_pos"] = participation_pos
        graph_prop_simps[st]["participation_neg"] = participation_neg
        #graph_prop_simps[st]["diversity_pos"] = diversity_pos # Over 64 nodes
        #graph_prop_simps[st]["diversity_neg"] = diversity_neg
        #pdb.set_trace()
        #graph_prop_simps[st]["median_node_strengths_pos"] = np.hstack(node_stren_pos)
        #graph_prop_simps[st]["median_node_strengths_neg"] = np.hstack(node_stren_neg)
        #graph_prop_simps[st]["total_node_strengths_pos"] = np.hstack(tot_node_stren_pos)
        #graph_prop_simps[st]["total_node_strengths_neg"] = np.hstack(tot_node_stren_neg)
        #graph_prop_simps[st]["gateway_coef_pos"] = gateway_coef_pos # Over 64 nodes
        #graph_prop_simps[st]["gateway_coef_neg"] = gateway_coef_neg
        graph_prop_simps[st]["module_degree_zscore"] = zscore
        #graph_prop_simps[st]["transitivity"] = transitivity

        print(len(norms),len(st_list_corr))

        nz_inds = np.unique(nz_inds)
        if len(norms) > len(st_list_corr):
            '''
            nz_inds = np.unique(nz_inds)
            diff_nz = np.diff(np.unique(nz_inds))
            ind1 = np.where(diff_nz>1)[0]
            if len(ind1)> 0:
                for ind2 in ind1:
                    ind_nz = ind1[0]+1
                    nz_inds[ind_nz:] = nz_inds[ind_nz:]-1
            '''
            graph_prop_simps[st]["st_list_corr_norm"] = np.array(norms)[nz_inds]
            graph_prop_simps[st]["total_amplitude"] = np.array(tot_amp)[nz_inds]
            graph_prop_simps[st]["average_amplitude"] = np.array(avg_amp)[nz_inds]
            graph_prop_simps[st]["percentage_active_sites"] = np.array(per_act_sit)[nz_inds]

            #graph_prop_simps[st]["clus_coef_pos"] = np.array(clus_coef_pos)[nz_inds]
            #graph_prop_simps[st]["clus_coef_neg"] = np.array(clus_coef_neg)[nz_inds]
            #graph_prop_simps[st]["loc_assort_pos"] = np.array(loc_ass_pos)[nz_inds]
            #graph_prop_simps[st]["loc_assort_neg"] = np.array(loc_ass_neg)[nz_inds]

            #graph_prop_simps[st]["names"] = np.array(names)[nz_inds]
        else:
            graph_prop_simps[st]["st_list_corr_norm"] = np.array(norms)
            #graph_prop_simps[st]["clus_coef_pos"] = np.array(clus_coef_pos)
            #graph_prop_simps[st]["clus_coef_neg"] = np.array(clus_coef_neg)
            graph_prop_simps[st]["total_amplitude"] = np.array(tot_amp)
            graph_prop_simps[st]["average_amplitude"] = np.array(avg_amp)
            graph_prop_simps[st]["percentage_active_sites"] = np.array(per_act_sit)

        if len(loc_ass_pos) > len(st_list_corr):
            graph_prop_simps[st]["local_assortativity_pos_whole"] = np.array(loc_ass_pos)[nz_inds]
            #graph_prop_simps[st]["loc_assort_neg"] = np.array(loc_ass_neg)[nz_inds]
        else:
            graph_prop_simps[st]["local_assortativity_pos_whole"] = np.array(loc_ass_pos)
            #graph_prop_simps[st]["loc_assort_neg"] = np.array(loc_ass_neg)


        if len(graph_properties[st]['names']) > len(st_list_corr):
            graph_prop_simps[st]["names"] = np.array(graph_properties[st]['names'])[nz_inds]
        else:
            graph_prop_simps[st]["names"] = np.array(graph_properties[st]['names'])

        if num_or_size == "num":
            ind_prop = 1
        elif num_or_size == "size":
            ind_prop = 2
        for k in np.arange(0,len(gammas)):

            if k == 0 or k == 7:
                if k == 0:
                    t1.plot(np.array(st_list_cov)[:,:,k][:,0],np.array(st_list_cov)[:,:,k][:,ind_prop],'.',color=cmaps[i](k+4),alpha=0.5,markersize=10)
                    t1.plot(np.median(np.array(st_list_cov)[:,:,k][:,0]),np.median(np.array(st_list_cov)[:,:,k][:,ind_prop]),'*',color=cmaps[i](k+2),markersize=20,label=st)
                    t2.plot(np.array(st_list_corr)[:,:,k][:,0],np.array(st_list_corr)[:,:,k][:,ind_prop],'.',color=cmaps[i](k+4),alpha=0.5,markersize=10)
                    t2.plot(np.median(np.array(st_list_corr)[:,:,k][:,0]),np.median(np.array(st_list_corr)[:,:,k][:,ind_prop]),'*',color=cmaps[i](k+4),markersize=20,label=st)
                else:
                    t1.plot(np.array(st_list_cov)[:,:,k][:,0],np.array(st_list_cov)[:,:,k][:,ind_prop],'.',color=cmaps[i](k+4),alpha=0.5,markersize=8)
                    t2.plot(np.array(st_list_corr)[:,:,k][:,0],np.array(st_list_corr)[:,:,k][:,ind_prop],'.',color=cmaps[i](k+4),alpha=0.5,markersize=8)
                    t1.plot(np.median(np.array(st_list_cov)[:,:,k][:,0]),np.median(np.array(st_list_cov)[:,:,k][:,ind_prop]),'*',color=cmaps[i](k+2),markersize=20)
                    #check_corr = anal.find_line_fit(np.array(st_list_corr)[:,:,k][:,0],np.array(st_list_corr)[:,:,k][:,1],1)
                    t2.plot(np.median(np.array(st_list_corr)[:,:,k][:,0]),np.median(np.array(st_list_corr)[:,:,k][:,ind_prop]),'*',color=cmaps[i](k+4),markersize=20)
            #check_cov = anal.find_line_fit(np.array(st_list_cov)[:,:,k][:,0],np.array(st_list_cov)[:,:,k][:,1],1)
            #t1.plot(np.array(st_list_cov)[:,:,k][:,0], check_cov,'-',color=cmaps[i](k+2))
            '''
            temp_dict["modularity_index"].append(np.array(st_list_corr)[:,:,k][:,0])
            temp_dict["largest_module_size"].append(np.array(st_list_corr)[:,:,k][:,2])
            temp_dict["num_modules"].append(np.array(st_list_corr)[:,:,k][:,1])
            temp_dict["num_modules_greater_than_"+str(percentile)].append(np.array(st_list_corr)[:,:,k][:,3])
            #temp_dict["within_module_distance"].append(np.array(graph_prop_simps[st]["within_module_dists"])[:,k])
            #med_within = [np.median(sl) for sl in np.array(graph_prop_simps[st]["within_module_dists"])[:,k]]
            med_within = [np.percentile(sl,50) for sl in np.array(graph_prop_simps[st]["within_module_dists"])[:,k]]
            temp_dict["median_within_module_distance"].append(med_within)
            #temp_dict["between_module_distance"].append(np.array(graph_prop_simps[st]["between_module_dists"])[:,k])
            #med_between = [np.median(sl) for sl in np.array(graph_prop_simps[st]["between_module_dists"])[:,k]]
            med_between = [np.percentile(sl,50) for sl in np.array(graph_prop_simps[st]["between_module_dists"])[:,k]]
            temp_dict["median_between_module_distance"].append(med_between)

            temp_dict["ratio_within_between"].append(np.array(med_within)/np.array(med_between))
            ''' 

            nz_inds = np.unique(nz_inds)
            temp_dict["gamma"].append([ np.round(gammas[k],2) for i2 in np.arange(0,len(np.array(st_list_corr)[:,:,k][:,0]))])
            if len(norms) > len(st_list_corr):
                '''
                nz_inds = np.unique(nz_inds)
                diff_nz = np.diff(np.unique(nz_inds))
                ind1 = np.where(diff_nz>1)[0]
                if len(ind1)> 0:
                    ind_nz = ind1[0]+1
                    nz_inds[ind_nz:] = nz_inds[ind_nz:]-1
                '''
                temp_dict["norms"].append(np.array(norms)[nz_inds])
                temp_dict["total_amplitude"].append(np.array(tot_amp)[nz_inds])
                temp_dict["average_amplitude"].append(np.array(avg_amp)[nz_inds])
                temp_dict["percentage_active_sites"].append(np.array(per_act_sit)[nz_inds])
                #temp_dict["clustering_coef_pos"].append(np.array(clus_coef_pos)[nz_inds])
                #temp_dict["clustering_coef_neg"].append(np.array(clus_coef_neg)[nz_inds])
                temp_dict["participation_pos"].append(np.array(graph_prop_simps[st]["participation_pos"])[nz_inds,k])
                temp_dict["participation_neg"].append(np.array(graph_prop_simps[st]["participation_neg"])[nz_inds,k])
                #temp_dict["diversity_pos"].append(np.array(graph_prop_simps[st]["diversity_pos"])[nz_inds,k])
                #temp_dict["diversity_neg"].append(np.array(graph_prop_simps[st]["diversity_neg"])[nz_inds,k])
                #temp_dict["gateway_coef_pos"].append(np.array(graph_prop_simps[st]["gateway_coef_pos"])[nz_inds,k])
                #temp_dict["gateway_coef_neg"].append(np.array(graph_prop_simps[st]["gateway_coef_neg"])[nz_inds,k])
                temp_dict["module_degree_zscore"].append(np.array(graph_prop_simps[st]["module_degree_zscore"])[nz_inds,k])
                #temp_dict["transitivity"].append(np.array(graph_prop_simps[st]["transitivity"])[nz_inds])

                #temp_dict["median_node_strengths_pos"].append(np.array(graph_prop_simps[st]["median_node_strengths_pos"])[nz_inds])
                #temp_dict["median_node_strengths_neg"].append(np.array(graph_prop_simps[st]["median_node_strengths_neg"])[nz_inds])
                #temp_dict["total_node_strengths_pos"].append(np.array(graph_prop_simps[st]["total_node_strengths_pos"])[nz_inds])
                #temp_dict["total_node_strengths_neg"].append(np.array(graph_prop_simps[st]["total_node_strengths_neg"])[nz_inds])

                #temp_dict["local_assortativity_pos"].append(np.array(graph_prop_simps[st]["loc_assort_pos"])[nz_inds])
                #temp_dict["local_assortativity_neg"].append(np.array(graph_prop_simps[st]["loc_assort_neg"])[nz_inds])

            else:
                temp_dict["norms"].append(np.array(norms))
                temp_dict["total_amplitude"].append(np.array(tot_amp))
                temp_dict["average_amplitude"].append(np.array(avg_amp))
                temp_dict["percentage_active_sites"].append(np.array(per_act_sit))

                #temp_dict["clustering_coef_pos"].append(clus_coef_pos)
                #temp_dict["clustering_coef_neg"].append(clus_coef_neg)
                temp_dict["participation_pos"].append(np.array(graph_prop_simps[st]["participation_pos"])[:,k])
                temp_dict["participation_neg"].append(np.array(graph_prop_simps[st]["participation_neg"])[:,k])
                #temp_dict["diversity_pos"].append(np.array(graph_prop_simps[st]["diversity_pos"])[:,k])
                #temp_dict["diversity_neg"].append(np.array(graph_prop_simps[st]["diversity_neg"])[:,k])
                #temp_dict["gateway_coef_pos"].append(np.array(graph_prop_simps[st]["gateway_coef_pos"])[:,k])
                #temp_dict["gateway_coef_neg"].append(np.array(graph_prop_simps[st]["gateway_coef_neg"])[:,k])
                temp_dict["module_degree_zscore"].append(np.array(graph_prop_simps[st]["module_degree_zscore"])[:,k])
                #temp_dict["transitivity"].append(np.array(graph_prop_simps[st]["transitivity"]))

                #temp_dict["median_node_strengths_pos"].append(np.array(graph_prop_simps[st]["median_node_strengths_pos"]))
                #temp_dict["median_node_strengths_neg"].append(np.array(graph_prop_simps[st]["median_node_strengths_neg"]))
                #temp_dict["total_node_strengths_pos"].append(np.array(graph_prop_simps[st]["total_node_strengths_pos"]))
                #temp_dict["total_node_strengths_neg"].append(np.array(graph_prop_simps[st]["total_node_strengths_neg"]))

            
            if len(names) > len(st_list_corr):
                temp_dict["names"].append(np.array(graph_prop_simps[st]["names"])[nz_inds])
            else:
                temp_dict["names"].append(np.array(graph_prop_simps[st]["names"]))
            
            temp_dict["local_assortativity_pos_whole"].append(np.array(graph_prop_simps[st]["local_assortativity_pos_whole"]))
            #temp_dict["local_assortativity_neg"].append(np.array(graph_prop_simps[st]["loc_assort_neg"]))

            count+=len(np.array(st_list_corr)[:,:,k][:,0]) 
        
        #if dat_type == "subtype":
        #temp_dict["subtype"].append( [st for i3 in np.arange(0,count)])
        temp_dict[dat_type].append( [st for i3 in np.arange(0,count)])

        #t1.set_xlim(0,9)
        #t2.set_xlim(0,9)

        print(st)
        print(len(st_list_cov))
        print(len(st_list_corr))


        #t1.plot(np.mean(np.array(st_list_cov)[:,0]),np.mean(np.array(st_list_cov)[:,1]),'*',label=st,color=cmap1(i),alpha=1.0,markersize=25)
        #t2.plot(np.mean(np.array(st_list_corr)[:,0]),np.mean(np.array(st_list_corr)[:,1]),'*',label=st,color=cmap1(i),alpha=1.0,markersize=25)
        #print(np.mean(np.array(st_list)[:,0]),np.mean(np.array(st_list)[:,1]))



    graph_prop_df["modularity_index"] = np.hstack(temp_dict["modularity_index"])
    #graph_prop_df["largest_module_size"] = np.hstack(temp_dict["largest_module_size"])
    #graph_prop_df["num_modules"] = np.hstack(temp_dict["num_modules"])
    #graph_prop_df["norms"] = np.hstack(temp_dict["norms"])
    graph_prop_df["total_amplitude"] = np.hstack(temp_dict["total_amplitude"])
    graph_prop_df["average_amplitude"] = np.hstack(temp_dict["average_amplitude"])
    graph_prop_df["percentage_active_sites"] = np.hstack(temp_dict["percentage_active_sites"])

    #graph_prop_df["num_modules_greater_than_"+str(percentile)] = np.hstack(temp_dict["num_modules_greater_than_"+str(percentile)])
    #graph_prop_df["within_module_distance"] = np.hstack(temp_dict["within_module_distance"])
    #graph_prop_df["median_within_module_distance"] = np.hstack(temp_dict["median_within_module_distance"])
    #graph_prop_df["between_module_distance"] = np.hstack(temp_dict["between_module_distance"])
    #graph_prop_df["median_between_module_distance"] = np.hstack(temp_dict["median_between_module_distance"])
    #graph_prop_df["ratio_within_between"] = np.hstack(temp_dict["ratio_within_between"])
    #graph_prop_df["clustering_coef_pos"] = np.hstack(temp_dict["clustering_coef_pos"])
    graph_prop_df["participation_pos"] = np.hstack(temp_dict["participation_pos"])
    #graph_prop_df["diversity_pos"] = np.hstack(temp_dict["diversity_pos"])
    #graph_prop_df["gateway_coef_pos"] = np.hstack(temp_dict["gateway_coef_pos"])
    graph_prop_df["local_assortativity_pos_whole"] = np.hstack(temp_dict["local_assortativity_pos_whole"])
    #graph_prop_df["clustering_coef_neg"] = np.hstack(temp_dict["clustering_coef_neg"])
    graph_prop_df["participation_neg"] = np.hstack(temp_dict["participation_neg"])
    #graph_prop_df["diversity_neg"] = np.hstack(temp_dict["diversity_neg"])
    #graph_prop_df["gateway_coef_neg"] = np.hstack(temp_dict["gateway_coef_neg"])
    #graph_prop_df["median_node_strengths_pos"] = np.hstack(temp_dict["median_node_strengths_pos"])
    #graph_prop_df["median_node_strengths_neg"] = np.hstack(temp_dict["median_node_strengths_neg"])
    #graph_prop_df["total_node_strengths_pos"] = np.hstack(temp_dict["total_node_strengths_pos"])
    #graph_prop_df["total_node_strengths_neg"] = np.hstack(temp_dict["total_node_strengths_neg"])
    graph_prop_df["module_degree_zscore"] = np.hstack(temp_dict["module_degree_zscore"])
    #graph_prop_df["transitivity"] = np.hstack(temp_dict["transitivity"])
    graph_prop_df["names"] = np.hstack(temp_dict["names"])
    #pdb.set_trace()


    #graph_prop_df["local_assortativity_neg"] = np.hstack(temp_dict["local_assortativity_neg"])
    '''
    graph_prop_df["ratio_cluster_coef_pos_neg"] = graph_prop_df["clustering_coef_pos"]/graph_prop_df["clustering_coef_neg"]
    graph_prop_df["ratio_participation_pos_neg"] = graph_prop_df["participation_pos"]/graph_prop_df["participation_neg"]
    graph_prop_df["ratio_diversity_pos_neg"] = graph_prop_df["diversity_pos"]/graph_prop_df["diversity_neg"]
    graph_prop_df["ratio_median_node_strengths_pos_neg"] = np.abs(graph_prop_df["median_node_strengths_pos"])/np.abs(graph_prop_df["median_node_strengths_neg"])
    graph_prop_df["ratio_total_node_strengths_pos_neg"] = np.abs(graph_prop_df["total_node_strengths_pos"])/np.abs(graph_prop_df["total_node_strengths_neg"])
    graph_prop_df["ratio_gateway_coef_pos_neg"] = graph_prop_df["gateway_coef_pos"]/graph_prop_df["gateway_coef_neg"]
    '''
    #graph_prop_df["ratio_local_assortativity_pos_neg"] = graph_prop_df["local_assortativity_pos"]/graph_prop_df["local_assortativity_neg"]



    graph_prop_df["gamma"] = np.hstack(temp_dict["gamma"])

    graph_prop_df[dat_type] = np.hstack(temp_dict[dat_type])
    #graph_prop_df["all_modules_sizes"] = np.hstack(temp_dict["all_modules_sizes"])

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


'''
t1.legend(loc='best',prop={'size':10,'weight':'bold'})
t2.legend(loc='best',prop={'size':10,'weight':'bold'})
t2.set_xlabel("Modularity index",fontsize=15,fontweight='bold')
if num_or_size == "num":
    t1.set_ylabel("Number of modules",fontsize=15,fontweight='bold')
    t2.set_ylabel("Number of modules",fontsize=15,fontweight='bold')
elif num_or_size == "size":
    t1.set_ylabel("Largest module size",fontsize=15,fontweight='bold')
    t2.set_ylabel("Largest module size",fontsize=15,fontweight='bold')

fig1.subplots_adjust(left=0.1,top = 0.96,bottom=0.1)
if data_type == "subtype":
    fig1.savefig(fig_target_dir+"graph_props_subtypes_"+num_or_size+".png")
elif data_type == "development":
    fig1.savefig(fig_target_dir+"graph_props_subtypes_days_"+num_or_size+".png")
'''

if data_type == "subtype":
    files1 = glob.glob(data_target_dir+"graph_properties_pandas_for_behav_[0-9]*.csv")
    files2 = glob.glob(data_target_dir+"graph_properties_pandas_[0-9]*.csv")
elif data_type == "development":
    files1 = glob.glob(data_target_dir+"graph_properties_pandas_for_behav_days_[0-9]*.csv")
    files2 = glob.glob(data_target_dir+"graph_properties_pandas_days_[0-9]*.csv")

def merge_df_seeds(files):
    pdb.set_trace()
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
    #merge_df1 = merge_df_seeds(files1)
    #merge_df1.to_csv(data_dir+"graph_properties_pandas_for_behav_all.csv") # everything

    merge_df2 = merge_df_seeds(files2)
    merge_df2.to_csv(data_dir+"graph_properties_pandas_all.csv") # nonan

elif data_type == "development":
    #merge_df1 = merge_df_seeds(files1)
    #merge_df1.to_csv(data_dir+"graph_properties_pandas_for_behav_days_all.csv")

    merge_df2 = merge_df_seeds(files2)
    merge_df2.to_csv(data_dir+"graph_properties_pandas_days_all.csv")




if data_type == "subtype":
    post_fix = ""
elif data_type == "development":
    post_fix = "_days_"
    #pdb.set_trace()




