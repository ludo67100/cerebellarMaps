import sys
import os
import numpy as np


cwd = os.getcwd()
os.chdir("../Behavior/")

sys.path.append(cwd)
cmd ='python glm_lcls.py auc_early'
os.system(cmd) 

