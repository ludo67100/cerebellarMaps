from itertools import combinations, product
import numpy as np
import re
import pickle

DATA_TARGET_DIR = "data/"
FIG_TARGET_DIR = "figs/"
SUBTYPES = ["EC","ENR1","ENR2","ES","LC","LS","WT"]
DEVELOPMENT =["P9P10","P12P13","P14P18","P30P40"]


Fig2_panel_name = dict({"modularity_index":"H","participation_pos":"I","module_degree_zscore":"J","local_assortativity_pos_whole":"K"})

#SEEDS = list(np.random.randint(1,99999999,20))
#pickle.dump(SEEDS,open(DATA_TARGET_DIR+"seeds.pickle","wb"))
SEEDS = list(pickle.load(open(DATA_TARGET_DIR+"seeds.pickle","rb")))[:5]
print(SEEDS)

rule all:
	input:
		DATA_TARGET_DIR+"meta_data.csv",
		DATA_TARGET_DIR+"data_2d_maps.pickle",
		expand(DATA_TARGET_DIR+"covariance_maps_norm_{seed}.pickle",seed=SEEDS),
		expand(FIG_TARGET_DIR+"corr_maps_rearranged_{st}_norm_{seed}.png",st=SUBTYPES,seed=SEEDS),
		DATA_TARGET_DIR+"meta_data_days.csv",
		DATA_TARGET_DIR+"data_2d_maps_days.pickle",	
		expand(DATA_TARGET_DIR+"covariance_maps_days_norm.pickle"),
		expand(DATA_TARGET_DIR+"graph_properties_days_norm_{seed}.pickle",seed=SEEDS),
		expand(FIG_TARGET_DIR+"corr_maps_rearranged_{st}_norm_{seed}.png",st=DEVELOPMENT,seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_days_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_days_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_days_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv",
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_all.csv",
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["modularity_index"])		






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


rule read_data_development:
	input:
		"For Paper/EPHYS/Development_Dataset/",
	output:
		DATA_TARGET_DIR+"meta_data_days.csv",
		DATA_TARGET_DIR+"data_2d_maps_days.pickle"	
	run:
		shell("python Development\ Dataset/Read_maps/read_data.py")
	
rule calc_graph_features_development:
	input:
		DATA_TARGET_DIR+"data_2d_maps_days.pickle",
		DATA_TARGET_DIR+"meta_data_days.csv"	
	output:
		expand(DATA_TARGET_DIR+"covariance_maps_days_norm.pickle"),
		expand(DATA_TARGET_DIR+"graph_properties_days_norm_{seed}.pickle",seed=SEEDS),
		expand(FIG_TARGET_DIR+"corr_maps_rearranged_{st}_norm_{seed}.png",st=DEVELOPMENT,seed=SEEDS)

	run:
		shell("python Development\ Dataset/Calculate_graph_features/plot_correlation_maps_and_calculate_graph_features.py")



rule collate_graph_features_development_dataframe:
	input:
		expand(DATA_TARGET_DIR+"graph_properties_days_norm_{seed}.pickle",seed=SEEDS)
	output:
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_days_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_days_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_days_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv"
	run:
		shell("python common/combine_graph_props_seeds_pandas.py development")


rule collate_graph_features_adaptive_dataframe:
	input:
		expand(DATA_TARGET_DIR+"graph_properties_norm_{seed}.pickle",seed=SEEDS)
	output:
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_all.csv"
	run:
		shell("python common/combine_graph_props_seeds_pandas.py subtype")


rule Figure2_panelH:
	input:
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv"
	output:
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["modularity_index"])
	run:
		shell("python Figure\ 2/Figure2_PanelH.py")
			

