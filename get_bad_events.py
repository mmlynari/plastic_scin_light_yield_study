import glob, os
import numpy as np
import pandas as pd
from multiprocessing import Pool
import matplotlib.pyplot as plt

# read in csv file
def csv2pd(csv_file):
    data = pd.read_csv(csv_file,delimiter=",",names = ["time", "amp_sipm","amp_trigger"],skiprows=21)
    data["time"] = data["time"].multiply(10**9) # convert time to nanoseconds
    if len( data["amp_trigger"].dropna() ) == 0 : # trigger included in data set ?
        data = data.drop(["amp_trigger"],axis=1)
    else:
        data["amp_trigger"] = data["amp_trigger"].multiply(-1) # invert waveform
        data["amp_trigger"] = data["amp_trigger"].add(3.5) # shift to positive values
    return data

def get_event_nr(filename):
	tree_name = filename.split("/")[-1]
	event_nr = 0
	test = "bla"
	event_nr_str = tree_name.replace("tek","").replace("CH2.csv","")
	# event_nr_str = tree_name.replace("tek","").replace("ALL.csv","")
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
	
nameSource = "13Nov2023_10kPoints_1"
# nameSource = "15Nov2023_day_10kpoints_2"
# nameSource = "15Nov2023_night_10kpoints_3"
# nameSource = "16Nov2023_night_10kpoints_4"
# nameSource = "17Nov2023_night_10kpoints_5"
# nameSource = "19Nov2023_day_10kpoints_6"

print("check %s"%nameSource)

# get list of filenames
file_list = glob.glob("/eos/user/t/tilepmt/SiPM/data/csv/data_10kpoints/"+nameSource+"/*.csv")
#print(file_list)

# df_list = []
bad_events = []
counter = 0
for name in file_list:
	try:
		event_nr = get_event_nr(name)
		# df_list.append(csv2pd(name))
		df = csv2pd(name)
		n_bins = len(df)
		if n_bins < 10000:
			print("not enough bins: n_bins = %d"%n_bins)
			print(event_nr)
			bad_events.append(event_nr)
			# print(df_list[0].shape)
		pass
	except pd.errors.ParserError: 
		print("pandas.errors.ParserError")
		print(name)
		print(event_nr)
		bad_events.append(event_nr)
		pass
	except UnicodeDecodeError:
		print("UnicodeDecodeError")
		print(name)
		bad_events.append(event_nr)
		print(event_nr)
		pass
	counter += 1

bad_events = np.array(bad_events)
print("list of bad event files:")
print(sorted(bad_events))
print("number of bad event files: %d"%len(bad_events))
