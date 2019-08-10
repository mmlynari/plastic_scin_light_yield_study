import glob, os, sys
import numpy as np
import ROOT as r

from multiprocessing import Pool
import pandas as pd
import matplotlib.pyplot as plt
import uproot

####################
## PE SLICES PROB ##
####################

def pe_slice_prob(df,pe_val):

	#__ DATAFRAME OF PE SLICES ___

	slice_names = []
	for i in pe_val:
		slice_names.append("N_pe{0}".format(i))	

	sc = 0 # slice counter
	d_pe_collection = {}
	# sum
	for k in pe_val:
		if k == 0:
			d_pe_collection[slice_names[sc]] = df.iloc[:,0][(df.iloc[:,0]<k+0.5) ]
		else:
			d_pe_collection[slice_names[sc]] = df.iloc[:,0][ (df.iloc[:,0]>k-0.5) & (df.iloc[:,0]<k+0.5) ]
		sc = sc+1

	df_pe_slices = pd.DataFrame.from_dict(d_pe_collection)

	slice_prob = np.zeros( (1,len(pe_val)) )
		# sum
	sc = 0		
	for i in range(0,len(pe_val)):
		slice_prob[0,i] = len( df_pe_slices[slice_names[sc]].dropna() ) / len(df)
		sc = sc+1

	df_slice_sum_prob = pd.DataFrame(data=slice_prob.T,columns=df.columns)
	
	return df_slice_sum_prob 

#################
## PE CUM PROB ##
#################

def pe_cum_prob(df_pe_slices_prob,pe_val):
	#__ DATAFRAME OF CUMMULATIVE PE PROBABILITIES ____
	slice_cprob = np.zeros( (1,len(pe_val)) )
	# sum
	sc = 0		
	for i in range(0,len(pe_val)):
		slice_cprob[0,i] = df_pe_slices_prob.iloc[i:len(pe_val),0].sum()
		sc = sc+1

	df_slice_cprob = pd.DataFrame(data=slice_cprob.T,columns=df_pe_slices_prob.columns)
	return df_slice_cprob

################
## INITIALIZE ##
################

directory = os.path.join( os.getcwd() )

# source directory
comb_root_dir = str(directory+"/comb_root")
if not os.path.isdir(comb_root_dir):
    os.system("mkdir -p %s"%comb_root_dir)


# ___ LIST ___

run_list = [
	"cosmics_1to6_015C_20ns_isCalib1_newBL",
	"cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL"
	# "comsics_runs_1to6_40ns"
]
var_list = [
	"dc_charge_alt",
	"dc_charge_alt"
]

run_nr_list = [
	[1,6],
	[7,50]
]
cut_list = [
	0,
	-0.001
]
clmn_list = [
	"015C",
	"1325CS"
]
label_list = [
	"Hamamatsu S12571-015C, int. window 20 ns",
	"Hamamatsu S13360-1325CS, int. window 100 ns"
	]
# label_list = [
# 	"SiPM: S12571-015C, int. win. 20 ns, runs {:d} to {:d}".format(run_nr_list[0][0],run_nr_list[0][1]),
# 	"SiPM: S13360-1325CS, int. win. 100 ns, runs {:d} to {:d}".format(run_nr_list[1][0],run_nr_list[1][1])
# ]

#############
## ANALYZE ##
#############

#__ loop over ALL RUNS ____

run_count = 0

pe_val = ([0,1,2,3,4,5,6,7,8,9,10])

# store DC spectra
d_data = {}

# store prob results
d_prob = {}
d_cprob = {}
v_mean = []

for runName in run_list:
	tree = uproot.open(comb_root_dir+"/{0}.root".format(runName))["ntuple"]

	# variable = "dc_charge_alt"
	variable = var_list[run_count]
	variable2 = "bl_value"
	variable3 = "run_nr"
	run_nr_l = run_nr_list[run_count][0]
	run_nr_u = run_nr_list[run_count][1]

	npe_lim = 15

	# get charge branch
	df_data = tree.pandas.df([variable,variable2,variable3])
	df_charge = pd.DataFrame( data=df_data.iloc[:,0][ (df_data.iloc[:,1]<cut_list[run_count]) & (df_data.iloc[:,2]<=run_nr_u) & (df_data.iloc[:,2]>=run_nr_l ) ].values,columns = [variable] )

	# probability of PE slices 
	df_pe_slices_prob = pe_slice_prob(df_charge,pe_val)	

	# print("PE SLICES PROB")
	# print(df_pe_slices_prob.head(10))
	# print(df_pe_slices_prob.sum(axis=0))

	# cumulative probability of PE slices 
	df_pe_slice_cprob = pe_cum_prob(df_pe_slices_prob,pe_val)

	# print("CUMM. PROB")	
	# print(df_pe_slice_cprob)

	df_charge = df_charge.replace([np.inf, -np.inf], np.nan)
	mean_charge = df_charge.mean()

	#___ SAVE RESULTS ____
	d_data["{0}".format(clmn_list[run_count])] = df_charge.iloc[:,0]
	d_prob["prob_{0}".format(clmn_list[run_count])] = df_pe_slices_prob.iloc[:,0]
	d_cprob["cprob_{0}".format(clmn_list[run_count])] = df_pe_slice_cprob.iloc[:,0]
	v_mean.append( mean_charge )

	run_count = run_count + 1
df_data = pd.DataFrame.from_dict(d_data)
df_pe_val = pd.DataFrame(data=pe_val,columns=["Thr"])
df_prob = pd.DataFrame.from_dict(d_prob)
df_prob = pd.concat(([df_pe_val,df_prob]),axis=1)
df_cprob = pd.DataFrame.from_dict(d_cprob)
df_cprob = pd.concat(([df_pe_val,df_cprob]),axis=1)
df_mean = pd.DataFrame(np.array(v_mean))


###########
## PLOTS ##
###########

nBins = 200
x_min = -2
sum_lim = 25
ch_lim = 8

pe_def = "charge"

fig0, ax0 = plt.subplots( nrows=1, ncols=2, figsize=(9,5) )
# fig0.suptitle("")

ax_ch = []
axx_ch = []
for i in range(0,2):
	ax_ch.append(ax0[i].twinx())
	axx_ch.append(ax0[i])


ch_prob_color = "black"
ch_hist_color = "dodgerblue"
sum_prob_color = "black"
sum_hist_color = "red"
y_labelsize = 9
x_labelsize = 9
grid_alpha = 0.15
leg_fsize = 9

ch_npe_lim = 1

for i in range(0,2):
		
	axx = axx_ch[i]
	ax = ax_ch[i]

	nBins = int(((df_data[clmn_list[i]].max() - df_data[clmn_list[i]].min()) *10).round(0))
	
	ax.hist(df_data[clmn_list[i]].dropna().values, bins=nBins, histtype="stepfilled",color=ch_hist_color,alpha=0.8, log=False,label="Dark count dist.\nmean = {:1.2f} $N_{{pe}}$\nentries: {:d} ".format(df_mean.iloc[i,0],len(df_data[clmn_list[i]].dropna().values)),density=True,zorder=2)

	ax.set_title("{0}".format(label_list[i]),fontsize=10)
	ax.legend(loc="center right",fontsize=leg_fsize)
	ax.set_xticks((np.arange(x_min, ch_lim+1,1)))
	ax.set_xlim(x_min,ch_lim)
	ax.set_yscale("log")
	ax.tick_params(axis='y', labelcolor=ch_hist_color,labelsize=y_labelsize)

	axx.plot(pe_val,df_prob.iloc[:,i+1].values,color=ch_prob_color,linestyle="-",linewidth=1,ms=3,marker="s",label="\nDark count prob.\n$P(N_{{pe}}={:d})$ = {:1.1f} %".format(ch_npe_lim,df_prob.iloc[ch_npe_lim,i+1]*100),zorder=1)
	axx.plot(pe_val,df_cprob.iloc[:,i+1].values,color=ch_prob_color,linestyle=":",linewidth=1,ms=3,marker="o",label=("Dark count cum. prob.\n"+r"$P^{{\ast}}(N_{{pe}}={:d})$ = {:1.1f} %".format(ch_npe_lim,df_cprob.iloc[ch_npe_lim,i+1]*100)),zorder=1)

	axx.legend(loc="upper right",fontsize=leg_fsize)
	axx.grid(True,"both","both",alpha=grid_alpha,color=ch_prob_color)
	axx.tick_params(axis='y', labelcolor=ch_prob_color,labelsize=y_labelsize)
	axx.tick_params(axis='x',labelsize=y_labelsize)
	axx.set_yscale("log")

	axx.set_xlabel("Threshold [$N_{pe}$]",fontsize = 9)
	ax.set_ylabel("Dark count distribution",fontsize = 9,color=ch_hist_color)
	axx.set_ylabel(r"Dark count probability $P^{{(\ast)}}(N_{pe})$",fontsize = 9,color=ch_prob_color,va="baseline")

plt.subplots_adjust(left=0.07, right=0.93, top=0.9, bottom=0.1, hspace=0.35, wspace=0.45)

plt.show()