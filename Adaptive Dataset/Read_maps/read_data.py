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

# Read the raw data from here
data_dir = "/home/bahuguna/Work/Isope_data/Isope_data_cerebellar_maps/"
# Store processed data here
data_target_dir = "../../data/"
fig_target_dir = "../../figs/"

#Name of the subfolder with adaptive data
electrophys = "ELECTROPHY"

zone_names = ["B_contra","AX_contra","Alat_contra","Amed_contra","Amed_ipsi","Alat_ipsi","AX_ipsi","B_ipsi"]

zone_lims = [(-233,-133),(-133,-108),(-108,-58),(-58,0),(0,50),(50,100),(100,125),(125,235)]

data = pd.DataFrame()

# All adaptive types = CT, EC, LC, LS, ES, S/L-TR
subtypes = os.listdir(data_dir+electrophys)

subtype_list = []
rat_num_list = []
cell_num_list = []

data_2d = dict()


for st in subtypes:
    all_cells = os.listdir(data_dir+electrophys+"/"+st)
    cell_no = []
    rat_no = []
    data_2d[st] = dict()

    for cell in all_cells:
        if '(' in cell:
            rn,cn = cell.split('_')[1].split('(')
            cn = cn.split(')')[0]
            rat_no.append(rn)
            cell_no.append(cn)
            subtype_list.append(st)
            if rn not in list(data_2d[st].keys()):
                data_2d[st][rn] = dict()
                data_2d[st][rn]["name"] = [] 

            if os.path.exists(data_dir+electrophys+"/"+st+"/"+cell+"/"+cell+"_Amp_zscore_2D_OK.csv"):
                if cn not in list(data_2d[st][rn].keys()):
                    data_2d[st][rn][cn] = dict()

                #synaptic currents in the map, not zscored
                data_2d[st][rn][cn]["map_nz"] = pd.read_csv(data_dir+electrophys+"/"+st+"/"+cell+"/"+cell+"_Amp_2D_OK.csv",header=None)
                data_2d[st][rn][cn]["map"] = pd.read_csv(data_dir+electrophys+"/"+st+"/"+cell+"/"+cell+"_Amp_zscore_2D_OK.csv",header=None)      # Assumption is if rat number is same, at least the cell number is different, that is no duplicate folders for the same rat num and cell number
		
                arr =  np.array(data_2d[st][rn][cn]["map"])
                pos_centered = pd.read_csv(data_dir+electrophys+"/"+st+"/"+cell+"/"+cell+"_Positions_cp_centered_OK.csv",header=None)
                ind_contra = np.where(pos_centered<0)[0]
                ind_ipsi = np.where(pos_centered>0)[0]
                data_2d[st][rn][cn]["ind_ipsi"] = ind_ipsi

                data_2d[st][rn][cn]["ind_contra"] = ind_contra

                data_2d[st][rn][cn]["ind_zones"] = [  np.where(np.logical_and(pos_centered>=x[0],pos_centered<x[1])==True)[0]  for x in zone_lims] 

                data_2d[st][rn][cn]["pos_centered"] = pos_centered

                data_2d[st][rn]["name"].append((rn,cn))

    rat_num_list.append(rat_no)    
    cell_num_list.append(cell_no)


data["subtypes"] = subtype_list
data["rat_num"] = np.hstack(rat_num_list)
data["cell_num"] = np.hstack(cell_num_list)

data.to_csv(data_target_dir+"meta_data.csv")
pickle.dump(data_2d,open(data_target_dir+"data_2d_maps.pickle","wb"))















