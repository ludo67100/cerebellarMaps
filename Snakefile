from itertools import combinations, product
import numpy as np
import re
import pickle

DATA_TARGET_DIR = "data/"
FIG_TARGET_DIR = "figs/"
SUBTYPES = ["EC","ENR1","ENR2","ES","LC","LS","WT"]
DEVELOPMENT =["P9P10","P12P13","P14P18","P30P40"]
OPTIONS = ["n","y"]
POSTFIX = ["","_sub_ipsi_contra"]


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
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["modularity_index"]),
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["participation_pos"]),
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["module_degree_zscore"]),
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["local_assortativity_pos_whole"]),

		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_days_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_days_all.csv",
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_all.csv",
		expand(FIG_TARGET_DIR+"tsne_all_subtypes_{sub_ipsi_contra}_{st}_seeds.png",sub_ipsi_contra='y',st="subtype"),
		expand(DATA_TARGET_DIR+"behavior_features_pandas{pf}.csv",pf=POSTFIX),
		expand(DATA_TARGET_DIR+"graph_properties_with_behavior_pandas{pf}_all.csv",pf=POSTFIX),
		FIG_TARGET_DIR+"Accuracy_comparison_y_subtype.png",
		FIG_TARGET_DIR+"Confusion_matrix_random_forest_y_subtype.png",
		#expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_{pf}all.csv",pf=POSTFIX),
		expand(DATA_TARGET_DIR+"graph_properties_behavior_enr{pf}_all.csv",pf=POSTFIX),

		DATA_TARGET_DIR+"Predicted_actual_scatter_points_slope_n.csv",
		FIG_TARGET_DIR+"/enr_slope/"+"Predicted_actual_scatter_jittered_slope_n.png",	
		DATA_TARGET_DIR+"Predicted_actual_scatter_points_total_distance_n.csv",
		FIG_TARGET_DIR+"/enr_total_distance/"+"Predicted_actual_scatter_jittered_total_distance_n.png"	




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
		expand(FIG_TARGET_DIR+"corr_maps_rearranged_{st}_norm_{seed}.png",st=SUBTYPES,seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_norm_{seed}.pickle",seed=SEEDS)
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



rule collate_graph_features_development_ipsi_contra_dataframe:
	input:
		expand(DATA_TARGET_DIR+"graph_properties_days_norm_{seed}.pickle",seed=SEEDS)
	output:
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_days_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_days_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_days_all.csv"
	run:
		shell("python common/combine_graph_props_ipsi_contra_seeds_pandas.py development")


rule collate_graph_features_subtype_ipsi_contra_dataframe:
	input:
		expand(DATA_TARGET_DIR+"graph_properties_norm_{seed}.pickle",seed=SEEDS)
	output:
		expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_{seed}.csv",seed=SEEDS),
		expand(DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_{seed}.csv",seed=SEEDS),
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_sub_contra_ipsi_all.csv"
	run:
		shell("python common/combine_graph_props_ipsi_contra_seeds_pandas.py subtype")



rule Figure2_panelH:
	input:
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv"
	output:
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["modularity_index"])
	run:
		shell("python Figure\ 2/Figure2_PanelH.py")
			

rule Figure2_panelI:
	input:
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv"
	output:
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["participation_pos"])
	run:
		shell("python Figure\ 2/Figure2_PanelI.py")


rule Figure2_panelJ:
	input:
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv"
	output:
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["module_degree_zscore"])
	run:
		shell("python Figure\ 2/Figure2_PanelJ.py")

rule Figure2_panelK:
	input:
		DATA_TARGET_DIR+"graph_properties_pandas_days_all.csv"
	output:
		expand(FIG_TARGET_DIR+"Figure2_Panel{N}_development.png",N=Fig2_panel_name["local_assortativity_pos_whole"])
	run:
		shell("python Figure\ 2/Figure2_PanelK.py")


rule read_data_behavior:
	input:
		"For Paper/BEHAVIOR/BALANCE/Catwalk_Norm_Profiles_Cuff_Sham_Ctrl.xlsx",
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_all.csv",
		DATA_TARGET_DIR+"graph_properties_pandas_for_behav_sub_contra_ipsi_all.csv"
	output:
		expand(DATA_TARGET_DIR+"behavior_features_pandas{pf}.csv",pf=POSTFIX),
		expand(DATA_TARGET_DIR+"graph_properties_with_behavior_pandas{pf}_all.csv",pf=POSTFIX)
	run:
		for op in OPTIONS:
			shell("python Behavior/read_data_behavior.py {sub_ipsi_contra}".format(sub_ipsi_contra=op))



rule Figure5_panelC:
	input:
		DATA_TARGET_DIR+"graph_properties_with_behavior_pandas_sub_ipsi_contra_all.csv"
	output:
		expand(FIG_TARGET_DIR+"tsne_all_subtypes_{sub_ipsi_contra}_{st}_seeds.png",sub_ipsi_contra='y',st="subtype")
	run:
		shell("python Figure\ 5/Figure5_PanelC.py y")


rule Figure5_panelD:
	input:
		DATA_TARGET_DIR+"graph_properties_with_behavior_pandas_all.csv",
		DATA_TARGET_DIR+"graph_properties_with_behavior_pandas_sub_ipsi_contra_all.csv"	
	output:
		FIG_TARGET_DIR+"Accuracy_comparison_y_subtype.png",
		FIG_TARGET_DIR+"Confusion_matrix_random_forest_y_subtype.png"
	run:
		shell("python Figure\ 5/Figure5_PanelD.py y subtype")

rule read_enrichment:
	input:
		"For Paper/BEHAVIOR/ENRICHMENT/Enrichment.xlsx"
	output:
		#expand(DATA_TARGET_DIR+"graph_properties_pandas_for_behav_{pf}all.csv",pf=POSTFIX),
		expand(DATA_TARGET_DIR+"graph_properties_behavior_enr{pf}_all.csv",pf=POSTFIX)
	run:
		for op in OPTIONS:
			shell("python Behavior/read_enrichment.py {sub_ipsi_contra}".format(sub_ipsi_contra=op))



rule Figure6_PanelA:
	input:
		DATA_TARGET_DIR+"graph_properties_behavior_enr_all.csv"
	output:
		DATA_TARGET_DIR+"Predicted_actual_scatter_points_slope_n.csv",
		FIG_TARGET_DIR+"/enr_slope/"+"Predicted_actual_scatter_jittered_slope_n.png"	
	run:
		shell("python Figure\ 6/Figure6_PanelA.py")


rule Figure6_PanelB:
	input:
		DATA_TARGET_DIR+"graph_properties_behavior_enr_all.csv"
	output:
		DATA_TARGET_DIR+"Predicted_actual_scatter_points_total_distance_n.csv",
		FIG_TARGET_DIR+"/enr_total_distance/"+"Predicted_actual_scatter_jittered_total_distance_n.png"	
	run:
		shell("python Figure\ 6/Figure6_PanelB.py")



