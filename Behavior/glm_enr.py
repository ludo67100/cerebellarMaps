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
import scipy.stats as sp_st
import sys
import seaborn as sns
from matplotlib.lines import Line2D
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels as sms
from adjustText import adjust_text
from sklearn.model_selection import cross_val_score, StratifiedKFold, GridSearchCV, StratifiedShuffleSplit, ShuffleSplit

sys.path.append("../common/")
import analyze as anal
# Raw data
data_dir = "For Paper/BEHAVIOR/ENRICHMENT/"

data_target_dir = "../data/"
#data_target_dir = "/home/bahuguna/Work/Isope_data/"
fig_target_dir = "../figs/"

zone_names = ["B_contra","AX_contra","Alat_contra","Amed_contra","Amed_ipsi","Alat_ipsi","AX_ipsi","B_ipsi"]

prop1 = sys.argv[1] # behavioral feature - slope, total_distance 
ipsi_contra = "n" # Using whole graph properties

nsplits = 150 # Number of bootstrapping iterations 
test_ratio = 0.25 # Cross validation - training set - 75%, testing set - 25%

# Make a folder for every behavioral feature
if os.path.exists(fig_target_dir+"enr_"+prop1) == False:
    os.mkdir(fig_target_dir+"enr_"+prop1)
fig_target_dir = fig_target_dir+"enr_"+prop1+"/"


seed = np.random.randint(0,999999,1)
print(seed)
np.random.seed(seed)

fixed_slopes = "y"
if ipsi_contra == "n":
    independent_props = ["modularity_index","participation_pos","module_degree_zscore","local_assortativity_pos_whole"]


#------------------------------------------------------------------------------------------------------------------------#
# Average over the gammas and animals. We should have only one point per animal

def avg_over_gammas_and_animals(field,data,non_string_columns,metric='mean'):
    temp_dat1_wgp = dict()
    for k in non_string_columns:
        temp_dat1_wgp[k] = []

    if metric == "mean":
        sub_dat1_wgp = anal.mean_over_gammas(field,data,temp_dat1_wgp,non_string_columns)
    elif metric == "median":
        sub_dat1_wgp = anal.median_over_gammas(field,data,temp_dat1_wgp,non_string_columns)
    for k in non_string_columns:
        sub_dat1_wgp[k] = sub_dat1_wgp[k].astype('float')

    return sub_dat1_wgp

#------------------------------------------------------------------------------------------------------------------------#

if ipsi_contra == "n":
    graph_prop_enr_behav_whole = pd.read_csv(data_target_dir+"graph_properties_behavior_enr_all.csv") 

dat_wgp = graph_prop_enr_behav_whole

# Drop nan rows
dat_wgp = dat_wgp.dropna()

# Shorten the names
dat_wgp["short-names"] = [ x.split('-')[0]+"-"+x.split('-')[1]  for x in dat_wgp["names"] ]
unnamed =[x for x in dat_wgp.keys() if "Unnamed" in x]
non_string_columns = list(set(dat_wgp.keys())-set(unnamed+['names', 'subtype',"short-names","seed","gamma","mouse","time"]))

# Median over gammas 
sub_dat_wgp = avg_over_gammas_and_animals(["names","gamma"],dat_wgp,non_string_columns,"median") # Median over gammas
sub_dat_wgp["short-names"] = [ x.split('-')[0]+"-"+x.split('-')[1]  for x in sub_dat_wgp["names"] ]

# Mean over the cells
sub_dat1_wgp = avg_over_gammas_and_animals(["names"],sub_dat_wgp,non_string_columns,"mean") # Mean over cells
sub_dat11_wgp = sub_dat1_wgp[independent_props+["total_days","names"]+[prop1]]



# Change name of columns, for zone wise, because "-" is considered as subtraction
change_dict = dict()
if ipsi_contra == "zone_wise":

    independent_props1 = []
    for keys in sub_dat11_wgp.keys():
        if "modularity_index" in keys :
            independent_props1.append(keys)
            continue
        elif "-" not in keys:
            continue
        change_dict[keys] = keys.split('-')[0]+"_"+keys.split('-')[1]
        independent_props1.append(keys.split('-')[0]+"_"+keys.split('-')[1])
    sub_dat11_wgp = sub_dat11_wgp.rename(columns=change_dict)


subtypes = np.zeros(len(sub_dat11_wgp["total_days"]))
# STR = -1, LTR = 1
ind_str = np.where(np.array(sub_dat11_wgp["total_days"])==7.0)[0]
ind_ltr = np.where(np.array(sub_dat11_wgp["total_days"])==19.0)[0]

subtypes[ind_str] = -1
subtypes[ind_ltr] = 1
sub_dat11_wgp["total_days_val"] = subtypes


# Generate the formula string : Dep_var ~ Indep_var1 + Indep_var2 + .... + Indep_var1 * str_bias + Indep_var1 * ltr_bias ....
formula_string = prop1 +"~"

# Stand alone terms - Indep_var1, Indep_var2 ...
for i,k in enumerate(independent_props):
    formula_string = formula_string +k+"+"

# Add interaction terms - Indep_var1 * str_bias, Indep_var2 * str_bias ...
for i,k in enumerate(independent_props):
    if i < len(independent_props)-1:
        formula_string = formula_string +k+'*total_days_val'+"+"
    else:
        formula_string = formula_string +k+'*total_days_val'

formula_string = formula_string+"-total_days_val" # So that this is not used as in dependent variable
print(formula_string)



sub_dat11_wgp = sub_dat11_wgp.dropna()
sub_dat11_wgp = sub_dat11_wgp.sample(frac=1).reset_index(drop=True)

# Assign the subtype on the basis of number of days trained
sub_dat11_wgp["subtype"] = [ "S-TR" if x == 7.0 else "L-TR"   for x in sub_dat11_wgp["total_days"]]

# Seperate X (independent) and Y(dependent) for the GLMs
X = sub_dat11_wgp[independent_props+["total_days","total_days_val","subtype","names"]]
Y = sub_dat11_wgp[prop1]

# Generate bootstrapping splits using sklearn.model_selection.ShuffleSplit  
seed = np.random.randint(0,999999999)
skf = ShuffleSplit(n_splits=nsplits, test_size=test_ratio,train_size=1-test_ratio, random_state=seed)

corr_predicted_from_all = []  # correlation between actual and predicted values for the behavioral features for every bootstrapping iteration
corr_predicted_from_all_shuff = [] # correlation between actual and predicted values for the behavioral features for every bootstrapping iteration when the data is shuffled

lingress_vals =[] # linear regress scipy.stats.linregress values for every bootstrapping iteration when data was pooled( S-TR + L-TR)
lingress_vals7 =[] # scipy.stats.linregress values for every bootstrapping iteration for only S-TR mice
lingress_vals19 =[] # scipy.stats.linregress values for every bootstrapping iteration for only L-TR mice

lingress_vals_shuff =[] # linear regress scipy.stats.linregress values for every bootstrapping iteration when data was pooled ( S-TR + L-TR) and shuffled
lingress_vals_shuff7 =[]
lingress_vals_shuff19 =[]

# To store the predicted values of the behavioral feature for every bootstrapping iteration
predicted_actual_scatter_plots = pd.DataFrame(columns=["Predicted","Actual","correlation","pval","classifier-type","subtype","names"])
temp_dict = dict()
for k in predicted_actual_scatter_plots.keys():
    temp_dict[k] = []


for train_index, test_index in skf.split(X):
    data_temp_train = pd.DataFrame()
    data_temp_test = pd.DataFrame()
    data_temp_train_shuff = pd.DataFrame()
    data_temp_test_shuff = pd.DataFrame()

    #print(len(train_index),len(test_index))
    data_temp_train = X.iloc[train_index]
    data_temp_test = X.iloc[test_index]
    data_temp_train[prop1] = Y.iloc[train_index]
    data_temp_test[prop1] = Y.iloc[test_index]

    data_temp_train_shuff = data_temp_train.copy()
    data_temp_test_shuff = data_temp_test.copy()

    # Shuffle the dependent variable for the training set
    data_temp_train_shuff[prop1] = np.random.permutation(data_temp_train[prop1].values)
    data_temp_test_shuff[prop1] = data_temp_test[prop1].values

    # GLMs	
    md_pre = smf.glm(formula_string,data_temp_train, groups=data_temp_train["total_days_val"],family=sm.families.Gaussian())
    
    # If there is an error in some bootstrapping iteration, drop the analysis from this iteration and continue. This happens ocassionally 
    try:
        mdf_pre = md_pre.fit(method='nm',maxiter=500)
    except np.linalg.LinAlgError as err:
        continue

    print("============================Actual=======================================================")
    print( mdf_pre.summary())


    # GLMs on the shuffled data
    md_pre_shuff = smf.glm(formula_string,data_temp_train_shuff, groups=data_temp_train_shuff["total_days_val"],family=sm.families.Gaussian())
    try:
        mdf_pre_shuff = md_pre_shuff.fit(method='nm',maxiter=500)
    except np.linalg.LinAlgError as err:
        continue
    
    print("============================Shuffled=======================================================")
    print( mdf_pre_shuff.summary())

    # Predict the values of the behavioral features on the basis of GLMs coefficients for this bootstrapping iteration 
    y_pred_all =  np.array(mdf_pre.predict(data_temp_test))
    # Predict the values separately for S-TR/L-TR
    y_pred_all7 = np.array(mdf_pre.predict(data_temp_test.loc[data_temp_test["subtype"]=="S-TR"]))
    y_pred_all19 = np.array(mdf_pre.predict(data_temp_test.loc[data_temp_test["subtype"]=="L-TR"]))
    
    # Predicting when GLM was trained for shuffled values
    y_pred_all_shuff =  np.array(mdf_pre_shuff.predict(data_temp_test_shuff))
    y_pred_all_shuff7 = np.array(mdf_pre_shuff.predict(data_temp_test_shuff.loc[data_temp_test_shuff["subtype"]=="S-TR"]))
    y_pred_all_shuff19 = np.array(mdf_pre_shuff.predict(data_temp_test_shuff.loc[data_temp_test_shuff["subtype"]=="L-TR"]))

    # Pearson correlation between actual and predicted value of the behavioral feature     
    (c1,p1) = sp_st.pearsonr(np.array(data_temp_test[prop1]),y_pred_all)
    corr_predicted_from_all.append((c1,p1))
    # chance correlations
    (c1_s,p1_s) = sp_st.pearsonr(np.array(data_temp_test_shuff[prop1]),y_pred_all_shuff)
    corr_predicted_from_all_shuff.append((c1_s,p1_s))
    
    # scipy.stats.linregress values for pooled data, only S-TR, only L-TR 
    slope,intercept, r_value, p_value, std_err = sp_st.linregress(np.array(data_temp_test[prop1]),y_pred_all)
    slope2, intercept2, r_value2, p_value2, std_err2 = sp_st.linregress(np.array(data_temp_test_shuff[prop1]),y_pred_all_shuff)
    lingress_vals.append((slope,intercept,r_value,p_value,std_err))
    lingress_vals_shuff.append((slope2,intercept2,r_value2,p_value2,std_err2))

    if len(y_pred_all7) > 0:
	    slope7,intercept7, r_value7, p_value7, std_err7 = sp_st.linregress(np.array(data_temp_test.loc[data_temp_test["subtype"]=="S-TR"][prop1]),y_pred_all7)

	    slope7_shuff, intercept7_shuff, r_value7_shuff, p_value7_shuff, std_err7_shuff = sp_st.linregress(np.array(data_temp_test_shuff.loc[data_temp_test_shuff["subtype"]=="S-TR"][prop1]),y_pred_all_shuff7)
	    lingress_vals7.append((slope7,intercept7,r_value7,p_value7,std_err7))

	    lingress_vals_shuff7.append((slope7_shuff,intercept7_shuff,r_value7_shuff,p_value7_shuff,std_err7_shuff))

    if len(y_pred_all19) > 0:
        slope19,intercept19, r_value19, p_value19, std_err19 = sp_st.linregress(np.array(data_temp_test.loc[data_temp_test["subtype"]=="L-TR"][prop1]),y_pred_all19)

        slope19_shuff, intercept19_shuff, r_value19_shuff, p_value19_shuff, std_err19_shuff = sp_st.linregress(np.array(data_temp_test_shuff.loc[data_temp_test_shuff["subtype"]=="L-TR"][prop1]),y_pred_all_shuff19)

        lingress_vals19.append((slope19,intercept19,r_value19,p_value19,std_err19))

        lingress_vals_shuff19.append((slope19_shuff,intercept19_shuff,r_value19_shuff,p_value19_shuff,std_err19_shuff))


    temp_dict["Predicted"].append(np.array(y_pred_all))
    temp_dict["Actual"].append(np.array(data_temp_test[prop1]))
    temp_dict["correlation"].append([c1 for i in np.arange(len(y_pred_all))])
    temp_dict["pval"].append([p1 for i in np.arange(len(y_pred_all))])
    temp_dict["classifier-type"].append(["Actual" for i in np.arange(len(y_pred_all))])
    temp_dict["subtype"].append(data_temp_test["subtype"])	
    temp_dict["names"].append(data_temp_test["names"])	
 
    temp_dict["Predicted"].append(np.array(y_pred_all_shuff))
    temp_dict["Actual"].append(np.array(data_temp_test_shuff[prop1]))
    temp_dict["correlation"].append([c1_s for i in np.arange(len(y_pred_all_shuff))])
    temp_dict["pval"].append([p1_s for i in np.arange(len(y_pred_all_shuff))])
    temp_dict["classifier-type"].append(["Shuffled" for i in np.arange(len(y_pred_all_shuff))])
    temp_dict["subtype"].append(data_temp_test_shuff["subtype"])	
    temp_dict["names"].append(data_temp_test_shuff["names"])	


if ipsi_contra == "n":
    tit1 = "Whole graph properties"


# Store in pandas dataframe
for k in predicted_actual_scatter_plots.keys():
    predicted_actual_scatter_plots[k] = np.hstack(temp_dict[k])

# Store it for future use
predicted_actual_scatter_plots.to_csv(data_target_dir+"Predicted_actual_scatter_points_"+prop1+"_"+ipsi_contra+".csv")


# R-value distribution for shuffled and actual GLM fits
all_rvals_df = pd.DataFrame(columns=["value_type","value","data_type","subtype"])
temp_dict_rval = dict()

for k in all_rvals_df.keys():
    temp_dict_rval[k] = []

temp_dict_rval["value"].append(np.array(lingress_vals)[:,2])
temp_dict_rval["value_type"].append(["r_value" for i in np.arange(len(np.array(lingress_vals)[:,2]))])
temp_dict_rval["data_type"].append(["Actual" for i in np.arange(len(np.array(lingress_vals)[:,2]))])
temp_dict_rval["subtype"].append(["Pooled" for i in np.arange(len(np.array(lingress_vals)[:,2]))])

temp_dict_rval["value"].append(np.array(lingress_vals7)[:,2])
temp_dict_rval["value_type"].append(["r_value" for i in np.arange(len(np.array(lingress_vals7)[:,2]))])
temp_dict_rval["data_type"].append(["Actual" for i in np.arange(len(np.array(lingress_vals7)[:,2]))])
temp_dict_rval["subtype"].append(["S-TR" for i in np.arange(len(np.array(lingress_vals7)[:,2]))])

# R-value distribution for the shuffled data
temp_dict_rval["value"].append(np.array(lingress_vals_shuff)[:,2])
temp_dict_rval["value_type"].append(["r_value" for i in np.arange(len(np.array(lingress_vals_shuff)[:,2]))])
temp_dict_rval["data_type"].append(["Shuffled" for i in np.arange(len(np.array(lingress_vals_shuff)[:,2]))])
temp_dict_rval["subtype"].append(["Pooled" for i in np.arange(len(np.array(lingress_vals_shuff)[:,2]))])


temp_dict_rval["value"].append(np.array(lingress_vals_shuff7)[:,2])
temp_dict_rval["value_type"].append(["r_value" for i in np.arange(len(np.array(lingress_vals_shuff7)[:,2]))])
temp_dict_rval["data_type"].append(["Shuffled" for i in np.arange(len(np.array(lingress_vals_shuff7)[:,2]))])
temp_dict_rval["subtype"].append(["S-TR" for i in np.arange(len(np.array(lingress_vals_shuff7)[:,2]))])


# Store in pandas data frame
for k in all_rvals_df.keys():
    all_rvals_df[k] = np.hstack(temp_dict_rval[k])



palette = [ "darkorange","dimgrey","firebrick","darkkhaki"]

g2 = sns.jointplot(x="Actual",y="Predicted",hue="classifier-type",data=predicted_actual_scatter_plots.loc[predicted_actual_scatter_plots["classifier-type"]=="Actual"],kind='scatter', joint_kws={'alpha': 0.0},palette=[palette[0]],legend=False,height=10)

xticks = g2.fig.axes[0].get_xticks()
xticklabels = [ str(x)  for x in xticks]
df_temp = predicted_actual_scatter_plots.loc[predicted_actual_scatter_plots["classifier-type"]=="Actual"]
animal_labels =[]
label_pos =[]

# Plot (Predicted,actual) scatter points with different colors for S-TR and L-TR animals
for i,xt in enumerate(np.unique(df_temp["Actual"])):
    val = xt
    temp1 = df_temp.loc[df_temp["Actual"]==val]
    an_name = np.unique(temp1["names"])[0]
	
    if np.unique(temp1["subtype"])=="S-TR": 
    	xm,ym = anal.plot_scatter_plot_errorbars(temp1,["Actual", "Predicted"],palette[2],'o',g2.fig.axes[0],return_pts=True) # Error bars for predicted around median
    	animal_labels.append(an_name)    
    	label_pos.append((xm,ym))

    elif np.unique(temp1["subtype"])=="L-TR":
        xm,ym = anal.plot_scatter_plot_errorbars(temp1,["Actual", "Predicted"],palette[3],'o',g2.fig.axes[0],return_pts=True)
        animal_labels.append(an_name)    
        label_pos.append((xm,ym))

g2.fig.axes[0].set_xticks(xticks)
g2.fig.axes[0].set_xticklabels(xticklabels)

dat_act = predicted_actual_scatter_plots.loc[predicted_actual_scatter_plots["classifier-type"]=="Actual"]
dat_shuff = predicted_actual_scatter_plots.loc[predicted_actual_scatter_plots["classifier-type"]=="Shuffled"]
slope, intercept, r_value, p_value, std_err = sp_st.linregress(dat_act["Actual"],dat_act["Predicted"])
slope2, intercept2, r_value2, p_value2, std_err2 = sp_st.linregress(dat_shuff["Actual"],dat_shuff["Predicted"])


sns.regplot(x="Actual",y="Predicted",data=dat_act,scatter=False,line_kws={'linewidth':4.5,'label':"Actual - (r="+str(np.round(r_value,2))+", p="+str(np.round(p_value,4))+")"},color=palette[0],ax=g2.fig.axes[0])
sns.regplot(x="Actual",y="Predicted",data=dat_shuff,scatter=False,line_kws={'linewidth':4.5,'label':"Shuffled - (r="+str(np.round(r_value2,2))+", p="+str(np.round(p_value2,4))+")"},color=palette[1],ax=g2.fig.axes[0])

sns.regplot(x="Actual",y="Predicted",data=dat_act.loc[dat_act["subtype"]=="S-TR"],scatter=False,line_kws={'linewidth':4.5,'linestyle':'dashed'},color=palette[2],ax=g2.fig.axes[0])

lines1 = [Line2D([0], [0], color=c, linewidth=4.5) for c in palette[2:]]
lines2 = [Line2D([0], [0], color=c, linewidth=4.5) for c in palette[:2]]

g2.fig.legend(lines1,("S-TR", "L-TR"),bbox_to_anchor=(0.7,0.8),prop={'weight':'bold'})
g2.fig.legend(lines2,("Actual", "Shuffled"),bbox_to_anchor=(0.85,0.8),prop={'weight':'bold'})

# Inset showing R-value distributions
left, bottom, width, height = [0.18, 0.6, 0.25, 0.25]
ax2 = g2.fig.add_axes([left, bottom, width, height])

all_rvals_df_sub = all_rvals_df.loc[(all_rvals_df["value_type"]=="r_value")]
sns.pointplot(x="data_type",y="value",hue="subtype",data=all_rvals_df_sub,dodge=0.2,capsize=0.2,palette=[palette[0],palette[2]],ax=ax2,linestyles=['-','--'])
for x in ax2.get_xticklabels():
    x.set_fontsize(12)
    x.set_fontweight('bold')

ax2.legend().set_title("")
ax2.set_xlabel("")
ax2.set_ylabel("")
inset_xlims = ax2.get_xlim()
inset_ylims = ax2.get_ylim()

ax2.hlines(xmin=inset_xlims[0],xmax=inset_xlims[1],y=0,color='k',linestyle='dashed')

ax2.set_title("r-value")

g2.ax_marg_x.remove()
g2.ax_marg_y.remove()

g2.fig.axes[0].set_ylabel(g2.fig.axes[0].get_ylabel()+" "+prop1,fontsize=15,fontweight='bold')
g2.fig.axes[0].set_xlabel(g2.fig.axes[0].get_xlabel()+" "+prop1,fontsize=15,fontweight='bold')

xmin,xmax = g2.fig.axes[0].get_xlim()
xbuff = 0.1*(xmax-xmin)
g2.fig.axes[0].set_xlim(xmin-xbuff,xmax+xbuff)


if prop1 == "total_distance":
	ymin,ymax = g2.fig.axes[0].get_ylim()
	ybuff = 0.05*(ymax-ymin)
	g2.fig.axes[0].set_ylim(ymin-ybuff,ymax+ybuff)
elif prop1 == "slope":
	g2.fig.axes[0].set_ylim(0,100)

if prop1 == "total_distance":
    if ipsi_contra == "n":
        g2.fig.axes[0].set_ylim(-5000,30000)
    else:
        g2.fig.axes[0].set_ylim(-90000,110000)

g2.fig.subplots_adjust(top=0.93,left=0.1)
g2.fig.savefig(fig_target_dir+"Predicted_actual_scatter_jittered_"+prop1+"_"+ipsi_contra+".png")
g2.fig.savefig(fig_target_dir+"Predicted_actual_scatter_jittered_"+prop1+"_"+ipsi_contra+".pdf")






















    



