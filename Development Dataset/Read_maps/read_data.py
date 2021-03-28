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
# Store the processed data here
data_target_dir = "../../data/"
# Target directory for the figures
fig_target_dir = "../../figs/"


# Name of the subfolder with developmental data
development = "DEVELOPMENT"

zone_names = ["B_contra","AX_contra","Alat_contra","Amed_contra","Amed_ipsi","Alat_ipsi","AX_ipsi","B_ipsi"]
zone_lims = [(-233,-133),(-133,-108),(-108,-58),(-58,0),(0,50),(50,100),(100,125),(125,235)]

'''
B_contra : -233 to -133
AX_contra : -133 to -108
Alat_contra : - 108 to -58
Amed_contra : -58 to 0

Amed ipsi :0 to 50
Alat_ipis = 50 to 100
AX_ipsi : 100 to 125
B_ipis : 125 to 285
'''

data = pd.DataFrame()

# All developmental stages - P9P10, P12P13, P14P18, P30P40
days = os.listdir(data_dir+development)

day_list = []
rat_num_list = []
cell_num_list = []

data_2d = dict()


for dy in days:
    all_cells = os.listdir(data_dir+development+"/"+dy)
    cell_no = []
    rat_no = []
    data_2d[dy] = dict()
    
    #Dictionary structure - Development stage - animal_number - cell num
    for cell in all_cells:
        if '(' in cell:
            rn,cn = cell.split('(')[0],cell.split('(')[1].split(')')[0]
            rat_no.append(rn)
            cell_no.append(cn)
            day_list.append(dy)
            if rn not in list(data_2d[dy].keys()):
                data_2d[dy][rn] = dict()
                data_2d[dy][rn]["name"] = [] 

            if os.path.exists(data_dir+development+"/"+dy+"/"+cell+"/"+cell+"_Amp_zscore_2D_OK.csv"):
                if cn not in list(data_2d[dy][rn].keys()):
                    data_2d[dy][rn][cn] = dict()

                # Read the zscored amplitude map
                data_2d[dy][rn][cn]["map"] = pd.read_csv(data_dir+development+"/"+dy+"/"+cell+"/"+cell+"_Amp_zscore_2D_OK.csv",header=None)      # Assumption is if rat number is same, at least the cell number is different, that is no duplicate folders for the same rat num and cell number
                data_2d[dy][rn][cn]["map_nz"] = pd.read_csv(data_dir+development+"/"+dy+"/"+cell+"/"+cell+"_Amp_2D_OK.csv",header=None)      # Assumption is if rat number is same, at least the cell number is different, that is no duplicate folders for the same rat num and cell number

                # Read the positions at which the amplitudes were measured
                pos_centered = pd.read_csv(data_dir+development+"/"+dy+"/"+cell+"/"+cell+"_Positions_cp_centered_OK.csv",header=None)

                # Mark the positions that fall on the contralateral side
                ind_contra = np.where(pos_centered<0)[0]
                # Mark the positions that fall on the ipsilateral side
                ind_ipsi = np.where(pos_centered>0)[0]
                data_2d[dy][rn][cn]["ind_ipsi"] = ind_ipsi
                data_2d[dy][rn][cn]["ind_contra"] = ind_contra

                # Zone wise indices
                data_2d[dy][rn][cn]["ind_zones"] = [  np.where(np.logical_and(pos_centered>=x[0],pos_centered<x[1])==True)[0]  for x in zone_lims] 

                data_2d[dy][rn][cn]["pos_centered"] = pos_centered
                data_2d[dy][rn]["name"].append((rn,cn)) # Rat number, cell number
    rat_num_list.append(rat_no)    
    cell_num_list.append(cell_no)


data["days"] = day_list
data["rat_num"] = np.hstack(rat_num_list)
data["cell_num"] = np.hstack(cell_num_list)

data.to_csv(data_target_dir+"meta_data_days.csv")
pickle.dump(data_2d,open(data_target_dir+"data_2d_maps_days.pickle","wb"))















