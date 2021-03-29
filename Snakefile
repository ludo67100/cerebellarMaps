from itertools import combinations, product
import numpy as np
import re
import pickle

DATA_TARGET_DIR = "data/"
FIG_TARGET_DIR = "figs/"
SUBTYPES = ["EC","ENR1","ENR2","ES","LC","LS","WT"]
#SEEDS = list(np.random.randint(1,99999999,20))
#pickle.dump(SEEDS,open(DATA_TARGET_DIR+"seeds.pickle","wb"))
SEEDS = list(pickle.load(open(DATA_TARGET_DIR+"seeds.pickle","rb")))[:5]

rule all:
	input:
		DATA_TARGET_DIR+"meta_data.csv",
		DATA_TARGET_DIR+"data_2d_maps.pickle",
		expand(DATA_TARGET_DIR+"covariance_maps_norm_{seed}.pickle",seed=SEEDS),
		expand(FIG_TARGET_DIR+"corr_maps_rearranged_{st}_norm_{seed}.png",st=SUBTYPES,seed=SEEDS)

		

rule read_data_adaptive:
	input:
		"For Paper/EPHYS/Adaptive_Dataset/",
	output:
		DATA_TARGET_DIR+"meta_data.csv",
		DATA_TARGET_DIR+"data_2d_maps.pickle",
	run:
		shell("python Adaptive\ Dataset/Read_maps/read_data.py")

rule calc_graph_features_adaptive:
	input:
		DATA_TARGET_DIR+"meta_data.csv",
		DATA_TARGET_DIR+"data_2d_maps.pickle",
	output:
		expand(DATA_TARGET_DIR+"covariance_maps_norm_{seed}.pickle",seed=SEEDS),
		expand(FIG_TARGET_DIR+"corr_maps_rearranged_{st}_norm_{seed}.png",st=SUBTYPES,seed=SEEDS)
	run:
		shell("python Adaptive\ Dataset/Calculate_graph_features/plot_correlation_maps_and_calculate_graph_features.py")


	

