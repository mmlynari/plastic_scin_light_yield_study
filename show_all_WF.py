import glob, os
import numpy as np
import pandas as pd
from multiprocessing import Pool
import matplotlib.pyplot as plt

# read in csv file
def csv2pd(csv_file):
	data = pd.read_csv(csv_file,delimiter=",",names = ["time", "amp_sipm","amp_trigger"],skiprows=21, usecols = ["time", "amp_sipm"])
	data["time"] = data["time"].multiply(10**9)
	data["amp_sipm"] = data["amp_sipm"].multiply(10**3) 
	# data = data.drop("amp_trigger",axis=1) # do not return trigger information
	return data

def get_event_nr(filename):
	tree_name = filename.split("/")[-1]
	event_nr = 0
	test = "bla"
	event_nr_str = tree_name.replace("tek","").replace("ALL.csv","")
	# print(event_nr_str)
	if event_nr_str == "0000":
		event_nr_str = event_nr_str.replace(event_nr_str[:3],"")
		event_nr = int(event_nr_str)
		print(event_nr)
	else:
		event_nr_str = event_nr_str.lstrip("0")
		try:
			event_nr_str = event_nr_str.lstrip("0")
			event_nr = int(event_nr_str)
			# print(event_nr)
		except ValueError:
			print("ValueError at event {0}".format(event_nr))
			pass
	return event_nr

# nameSource = "cosmics_bottom_1"
# nameSource = "cosmics_bottom_2"
# nameSource = "cosmics_bottom_3"
# nameSource = "cosmics_middle_4"
# nameSource = "cosmics_middle_7"
# nameSource = "cosmics_middle_10"
nameSource = "cosmics_middle_15"

wave_dir = "./root/"+nameSource+"/waveforms"

# number of parallel threads, depends on CPU cores
pool = Pool(processes=8) 

# get list of filenames
file_list = glob.glob("./csv/"+nameSource+"/*.csv")
# print(file_list)

# go through list in parallel
df_list = pool.map(csv2pd, file_list)

# combine dataframes
combined_df = pd.concat(df_list, ignore_index=True)
combined_df = combined_df.replace([np.inf, -np.inf], np.nan)
# print(combined_df.shape )

fig0, axs0 = plt.subplots( ncols=1, figsize=(10,7) )
fig0.suptitle("waveform overlay\nrun: {0}".format(nameSource))
ax = axs0
hb = ax.hexbin(combined_df["time"].values, combined_df["amp_sipm"].values, gridsize=(200,80), bins="log", cmap="inferno")
ax.ticklabel_format(axis="x",style="scientific")

ax.set_xlabel("time [ns]")
ax.set_ylabel("amplitude [mV]")

cb = fig0.colorbar(hb, ax=ax)
cb.set_label('log10(#Entries)')

plt.subplots_adjust(left=0.1, right=1.0, top=0.9, bottom=0.1, hspace=0.35, wspace=0.45)

# plt.show()	

fig0.savefig(wave_dir+"/wf_overlay.pdf")

