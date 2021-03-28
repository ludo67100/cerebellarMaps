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
import bct
from collections import Counter 
import seaborn as sns
import scipy.spatial.distance as sp_sp_dist
from itertools import product
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import scipy.cluster.hierarchy as shc
import sklearn as skl
from sklearn.cluster import KMeans
gammas = np.round(np.arange(0.0,1.5,0.17),2)

def calc_modularity(mat,weighted=True,undirected=True): 
    gammas = np.arange(0.0,1.5,0.17)
    num_mods_list = []
    modularity_index_list = []
    ci_list = []  # community vector list
    for g in gammas:

        if undirected == True:
            mod = bct.modularity_louvain_und_sign(mat,gamma=g)
        else:
            mod = bct.modularity.modularity_louvain_dir(mat,gamma=g)
        num_mods = [Counter(mod[0])[x]  for x in Counter(mod[0]).keys() if Counter(mod[0])[x] > 1 ]
        ind_mods = [np.where(mod[0]==x)[0]  for x in Counter(mod[0]).keys() if Counter(mod[0])[x] > 1 ]
        modularity_index = mod[1]
        num_mods_list.append(num_mods)
        modularity_index_list.append(modularity_index)
        ci_list.append(mod[0])
    return gammas, num_mods_list, modularity_index_list,ci_list


def calc_local_assortativity_sign(mat):
    loc_pos, loc_neg = bct.local_assortativity_wu_sign(mat)

    return loc_pos, loc_neg


def calc_module_degree_zscore(mat,ci,undirected=True,median=True):
    if undirected == True:
        zscore = bct.centrality.module_degree_zscore(mat,ci,0) # undirected
    else:
        zscore = bct.centrality.module_degree_zscore(mat,ci,3) # directed graph in and out degree

    if median == True:
        return np.median(zscore)
    else:
        return zscore


# Participation coefficient is a measure of diversity of intermodular connections of individual nodes
# Ppos (Nx1 numpy.ndarray) – participation coefficient from positive weights
# Pneg (Nx1 numpy.ndarray) – participation coefficient from negative weights
def calc_participation_coef_sign(mat,ci_list,median=True,undirected=True):
   
    if undirected == True:
        med_participation_pos = []
        med_participation_neg = []

        for ci in ci_list:
            part = bct.participation_coef_sign(mat,ci)
            if median == True:
                med_participation_pos.append(np.median(part[0]))
                med_participation_neg.append(np.median(part[1]))
            else:
                med_participation_pos.append(part[0])
                med_participation_neg.append(part[1])

        return med_participation_pos, med_participation_neg
    else:
        part = bct.centrality.participation_coef(mat,ci_list,degree='out')

        return part
    

def get_re_arranged_matrix(label_comms,orig_mat):
    re_arranged_mat = np.copy(orig_mat)
    idx = np.argsort(label_comms)
    re_arranged_mat = re_arranged_mat[idx,:]
    re_arranged_mat = re_arranged_mat[:,idx]

    return re_arranged_mat


