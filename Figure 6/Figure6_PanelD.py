import sys
import os
import numpy as np


cwd = os.getcwd()
os.chdir("../Behavior/")

sys.path.append(cwd)
cmd ='python glm_lcls.py post_op15'
os.system(cmd) 

