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


# Raw data
data_dir = "/home/bahuguna/Work/Isope_data/Isope_data_cerebellar_maps/"
# 
data_target_dir = "../data/"
fig_target_dir = "../figs/"

electrophys = "ELECTROPHY"
behavior = "BEHAVIOR"

sub_ipsi_contra = sys.argv[1]

if sub_ipsi_contra == "n":
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_pandas_for_behav_all.csv")
else:
    graph_prop_df = pd.read_csv(data_target_dir+"graph_properties_pandas_for_behav_sub_contra_ipsi_all.csv")


behavior_features = pd.DataFrame(columns=["baseline","auc_early","auc_late","auc_global","post_op2","post_op2_rel","post_op33", "post_op33_rel", "tot_auc_pos", "tot_auc_neg","ratio_auc_pos_neg","num_switches","time_to_peak_wrt_post_op2","time_to_peak_wrt_baseline", "pos_neg_switch_slope","neg_pos_switch_slope","post_op14","post_op15","variance","names","subtype"])
temp_behav = dict()

y_features = list(behavior_features.keys())

for k in y_features:
    temp_behav[k] = []
                                                                                                                 
behavior_catwalk = pd.read_excel(data_dir+behavior+"/"+"Catwalk.xlsx")

animals_names = [x.split('_')[1] for x in list(behavior_catwalk["mouse"]) ]
temp_subtypes = [x.split('_')[0] for x in list(behavior_catwalk["mouse"]) ]

gammas = np.round(np.arange(0.0,1.5,0.17),2)

baseline_all = np.zeros((len(graph_prop_df),1))
auc_early_all = np.zeros((len(graph_prop_df),1))
auc_late_all = np.zeros((len(graph_prop_df),1))
auc_global_all = np.zeros((len(graph_prop_df),1))
post_op2_all = np.zeros((len(graph_prop_df),1))
post_op14_all = np.zeros((len(graph_prop_df),1))
post_op15_all = np.zeros((len(graph_prop_df),1))
post_op2_rel_all = np.zeros((len(graph_prop_df),1)) # Relative to baseline
post_op33_all = np.zeros((len(graph_prop_df),1))
post_op33_rel_all = np.zeros((len(graph_prop_df),1)) # Relative to baseline
tot_pos_auc_all = np.zeros((len(graph_prop_df),1))
tot_neg_auc_all = np.zeros((len(graph_prop_df),1))
ratio_pos_neg_auc_all = np.zeros((len(graph_prop_df),1))
num_switches_all = np.zeros((len(graph_prop_df),1))
time_peak_all = np.zeros((len(graph_prop_df),1))
time_peak_wrt_baseline_all = np.zeros((len(graph_prop_df),1))
pos_neg_switch_slope_all = np.zeros((len(graph_prop_df),1))
neg_pos_switch_slope_all = np.zeros((len(graph_prop_df),1))
ff_all = np.zeros((len(graph_prop_df),1))

graph_names_short = [ x.split('-')[0] if len(x.split('-')) <= 2 else x.split('-')[0]+"-"+x.split('-')[1] for x in graph_prop_df["names"]]

for i,an in enumerate(animals_names):
    catwalk_mouse = behavior_catwalk[behavior_catwalk["mouse"].str.contains(an)]
    baseline = float(catwalk_mouse["baseline"])
    auc_early = np.sum(np.array((catwalk_mouse[catwalk_mouse.keys()[2:4]]))) # algebraic value with sign
    auc_late = np.sum(np.array((catwalk_mouse[catwalk_mouse.keys()[4:8]])))
    auc_global = np.sum(np.array(np.abs(catwalk_mouse[catwalk_mouse.keys()[2:8]]))) # Total plasticity irrespective of direction
    post_op2 = np.unique(catwalk_mouse["post_op_2"]) # Absolute postop2
    post_op2_rel = (np.unique(catwalk_mouse["post_op_2"]) - baseline)/baseline  # Relative to baseline - even after zscoring does not give a good clustering 
    post_op33 = np.unique(catwalk_mouse["post_op_33"])
    post_op33_rel = (np.unique(catwalk_mouse["post_op_33"]) - baseline)/baseline
    post_op14 = np.unique(catwalk_mouse["post_op_14"]) # Absolute postop14
    post_op15 = np.unique(catwalk_mouse["post_op_15"]) # Absolute postop14

    ff = np.var(np.array(catwalk_mouse[catwalk_mouse.keys()[2:8]]))

    tot_pos_auc = np.sum(np.array(catwalk_mouse[catwalk_mouse.keys()[1:-1]])[np.array(catwalk_mouse[catwalk_mouse.keys()[1:-1]])>0])
    tot_neg_auc = np.sum(np.array(catwalk_mouse[catwalk_mouse.keys()[1:-1]])[np.array(catwalk_mouse[catwalk_mouse.keys()[1:-1]])<0])
    ratio_pos_neg_auc = np.abs(tot_pos_auc/tot_neg_auc)
    traj = np.sign(np.array(catwalk_mouse)[0][2:-1]) 
    num_switches = ((traj[:-1]*traj[1:])<0).sum()
    pos_points = np.array(catwalk_mouse[catwalk_mouse.keys()[1:-1]])[np.array(catwalk_mouse[catwalk_mouse.keys()[1:-1]])>0]
    
    traj_org = np.array(catwalk_mouse[catwalk_mouse.keys()[3:-1]])
    if len(np.where(traj_org > 0)[0]) > 0:
        time_to_peak = np.argmax(traj_org) # After pos_op_2
    else:
        time_to_peak = 0

    traj_orig_baseline = np.array(catwalk_mouse[catwalk_mouse.keys()[2:-1]])
    if len(np.where(traj_orig_baseline>0)[0]) > 0:
        time_to_peak_baseline = np.argmax(traj_orig_baseline) # After pos_op_2
    else:
        time_to_peak_baseline = 0
    if len(pos_points) > 0:
        peak_pos = np.max(pos_points)
    else:
        peak_pos = 0.0
    
    ind_switches = np.where(traj[:-1]*traj[1:]<0)
    slopes = np.array([ np.array(catwalk_mouse)[0][2:][i+1] - np.array(catwalk_mouse)[0][2:][i]    for i in np.where(traj[:-1]*traj[1:]<0)[0] if i < len(np.array(catwalk_mouse)[0][1:])])
    neg_pos_slope = [ np.max(slopes[slopes>0]) if len(slopes[slopes>0]) > 0 else 0][0]
    pos_neg_slope = [ np.min(slopes[slopes<0]) if len(slopes[slopes<0]) > 0 else 0][0]

     
    ind = np.where(np.array(graph_names_short)==an)[0]
    baseline_all[ind] = baseline
    auc_early_all[ind] = auc_early
    auc_late_all[ind] = auc_late
    auc_global_all[ind] = auc_global
    post_op2_all[ind] = post_op2
    post_op2_rel_all[ind] = post_op2_rel
    post_op33_all[ind] = post_op33
    post_op33_rel_all[ind] = post_op33_rel
    tot_pos_auc_all[ind] = tot_pos_auc
    tot_neg_auc_all[ind] = tot_neg_auc
    ratio_pos_neg_auc_all[ind] = ratio_pos_neg_auc
    num_switches_all[ind] = num_switches
    time_peak_all[ind] = time_to_peak
    time_peak_wrt_baseline_all[ind] = time_to_peak_baseline

    neg_pos_switch_slope_all[ind] = neg_pos_slope
    pos_neg_switch_slope_all[ind] = pos_neg_slope
    post_op14_all[ind] = post_op14
    post_op15_all[ind] = post_op15
    ff_all[ind] = ff

    temp_behav["names"].append(an)
    temp_behav["baseline"].append(baseline)
    temp_behav["auc_early"].append(auc_early)
    temp_behav["auc_late"].append(auc_late)
    temp_behav["auc_global"].append(auc_global)
    temp_behav["post_op2"].append(post_op2[0])
    temp_behav["post_op14"].append(post_op14[0])
    temp_behav["post_op15"].append(post_op15[0])
    temp_behav["post_op2_rel"].append(post_op2_rel[0])
    temp_behav["post_op33"].append(post_op33[0])
    temp_behav["post_op33_rel"].append(post_op33_rel[0])
    temp_behav["tot_auc_pos"].append(tot_pos_auc)
    temp_behav["tot_auc_neg"].append(tot_neg_auc)
    temp_behav["time_to_peak_wrt_post_op2"].append(time_to_peak)
    temp_behav["time_to_peak_wrt_baseline"].append(time_to_peak_baseline)
    temp_behav["ratio_auc_pos_neg"].append(ratio_pos_neg_auc)
    temp_behav["neg_pos_switch_slope"].append(neg_pos_slope)
    temp_behav["pos_neg_switch_slope"].append(pos_neg_slope)
    temp_behav["num_switches"].append(num_switches)
    temp_behav["variance"].append(ff)
    temp_behav["subtype"].append(temp_subtypes[i])





graph_prop_df["baseline"] = baseline_all
graph_prop_df["auc_early"] = auc_early_all
graph_prop_df["auc_late"] = auc_late_all
graph_prop_df["auc_global"] = auc_global_all
graph_prop_df["post_op2"] = post_op2_all
graph_prop_df["post_op14"] = post_op14_all
graph_prop_df["post_op15"] = post_op15_all
graph_prop_df["post_op2_rel"] = post_op2_rel_all
graph_prop_df["post_op33"] = post_op33_all
graph_prop_df["post_op33_rel"] = post_op33_rel_all
graph_prop_df["tot_auc_pos"] = tot_pos_auc_all
graph_prop_df["tot_auc_neg"] = tot_neg_auc_all
graph_prop_df["ratio_auc_pos_neg"] = ratio_pos_neg_auc_all
graph_prop_df["time_to_peak_wrt_post_op2"] = time_peak_all
graph_prop_df["time_to_peak_wrt_baseline"] = time_peak_wrt_baseline_all
graph_prop_df["neg_pos_switch_slope"] = neg_pos_switch_slope_all
graph_prop_df["pos_neg_switch_slope"] = pos_neg_switch_slope_all

graph_prop_df["num_switches"] = num_switches_all
graph_prop_df["variance"] = ff_all

for k in list(behavior_features.keys()):
    behavior_features[k] = temp_behav[k]

if sub_ipsi_contra == "n":
    behavior_features.to_csv(data_target_dir+"behavior_features_pandas.csv")
else:
    behavior_features.to_csv(data_target_dir+"behavior_features_pandas_sub_ipsi_contra.csv")

if sub_ipsi_contra == "n": 
    graph_prop_df.to_csv(data_target_dir+"graph_properties_with_behavior_pandas_all.csv")
else:
    graph_prop_df.to_csv(data_target_dir+"graph_properties_with_behavior_pandas_sub_ipsi_contra_all.csv")





