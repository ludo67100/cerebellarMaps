import sys
import os
import numpy as np


cwd = os.getcwd()
os.chdir("common/")

#sys.path.append("../Development\ Dataset/Plot_graph_features/")
sys.path.append(cwd)
cmd ='python plot_graph_props_final_figs.py mean local_assortativity_pos_whole development'
os.system(cmd) 

