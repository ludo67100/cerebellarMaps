from itertools import combinations, product
import numpy as np

DATA_TARGET_DIR = "data/"
FIG_TARGET_DIR = "figs/"


rule all:
	input:
		DATA_TARGET_DIR+"meta_data.csv",
		DATA_TARGET_DIR+"data_2d_maps.pickle",
	
		

rule read_data_adaptive:
	input:
		"For Paper/EPHYS/Adaptive_Dataset/",
	output:
		DATA_TARGET_DIR+"meta_data.csv",
		DATA_TARGET_DIR+"data_2d_maps.pickle",
	run:
		shell("python Adaptive\ Dataset/Read_maps/read_data.py")

		

