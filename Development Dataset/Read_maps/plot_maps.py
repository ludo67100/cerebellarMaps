import os
import glob
import numpy as np
import pylab as pl
import scipy.io as sio
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


# All developmental stages - P9P10, P12P13, P14P18, P30P40
days = os.listdir(data_dir+development)
data = pd.read_csv(data_target_dir+"meta_data_days.csv")
data_2d = pickle.load(open(data_target_dir+"data_2d_maps_days.pickle","rb"))

# Plot all the 2D maps in a subtype to get an overview
for dy in days:
    data_slice = data.loc[data["days"]==dy]
    num_subfigs = len(data_slice)
    fig = pl.figure(figsize=(12,16))
    rows = int(num_subfigs/2)
    if rows*2 < num_subfigs:
        rows = rows+1
    subfig_hands = []

    fig.suptitle("Days:"+dy,fontsize=15,fontweight='bold') 
    for i,(rn,cn) in enumerate(zip(data_slice["rat_num"],data_slice["cell_num"])):
        subfig_hands.append(fig.add_subplot(rows,2,i+1))
        if str(cn) in list(data_2d[dy][str(rn)].keys()):
            map_2d = data_2d[dy][str(rn)][str(cn)]["map"]
            r,c = np.shape(map_2d)
            subfig_hands[-1].pcolor(map_2d,cmap=cm.hot)

            for i1,iz in enumerate(data_2d[dy][str(rn)][str(cn)]["ind_zones"]):
                if len(iz) == 0:
                    continue
    
                subfig_hands[-1].vlines(x=np.max(iz),ymin=0,ymax=r,color='white',linestyles='dashed',linewidth=1.5)                   
                subfig_hands[-1].text(np.min(iz)+(np.max(iz)-np.min(iz))/2.5,r*0.2,zone_names[i1],color='magenta',rotation=90)


            #subfig_hands[-1].set_aspect(5)
        subfig_hands[-1].set_title("rat num:"+str(rn)+",cell num:"+str(cn),fontsize=12,fontweight='bold')
        if i < (num_subfigs-2):
            subfig_hands[-1].set_xticklabels([])



    fig.subplots_adjust(left = 0.05,right=0.96,wspace=0.2,hspace=0.35,bottom=0.06,top=0.95)
    fig.savefig(fig_target_dir+"maps_"+dy+".png")



