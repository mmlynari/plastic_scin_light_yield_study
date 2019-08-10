#!/bin/bash

work_data()
{
	arguments=$0
	python3 convert_hist_m_thread.py $arguments
}

# analysis parameters
setCalib=0
setEvents=10000
setWindow=100

data_list=(
	# S12571-015C
	# "--nanoSec $setWindow --isCalib $setCalib --events 2001 --name cosmics_bottom_1"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 2981 --name cosmics_bottom_2"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 10000 --name cosmics_bottom_3"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 10000 --name cosmics_middle_4"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 5545 --name cosmics_middle_5"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 10000 --name cosmics_middle_6"
	# "--nanoSec $setWindow --isCalib $setCalib --events 9937 --name muons_tile10_2mFibre_103" 
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 10000 --name muons_tile10_2mFibre_104" 
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 7975 --name muons_tile10_2mFibre_105" 
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9999 --name muons_tile10_05mFibre_107" 
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9631 --name muons_tile10_05mFibre_108" 
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 8520 --name muons_tile10_05mFibre_109" 


	

	# S13360-1325CS
	# "--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_7"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_8"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_9"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_10"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 5336 --name cosmics_middle_11"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_12"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_13"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --name cosmics_middle_14"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9972 --name cosmics_middle_15"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9932 --name cosmics_middle_16"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9932 --name cosmics_middle_16"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 7073 --name cosmics_middle_32"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 6429 --name cosmics_middle_33"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 3791 --name cosmics_middle_34"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 8833 --name cosmics_middle_35"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 5235 --name cosmics_middle_36"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 8485 --name cosmics_middle_37"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 7386 --name cosmics_middle_38"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9206  --name cosmics_middle_39"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 5115  --name cosmics_middle_40"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 6839 --name cosmics_middle_41"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 8286 --name cosmics_middle_42"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 3808 --name cosmics_middle_43"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 6821 --name cosmics_middle_44"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 4917 --name cosmics_middle_45"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9825 --name cosmics_middle_46"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9791 --name cosmics_middle_47"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 5716 --name cosmics_middle_48"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents  --name cosmics_middle_49"
	# "\n--nanoSec $setWindow --isCalib $setCalib --events 9242 --name cosmics_middle_50"
	"\n--nanoSec $setWindow --isCalib $setCalib --events 100 --name cosmics_middle_50"

	
	) 

export -f work_data

echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "work_data"




#################
### DEBUGGING ###
#################

# Test functions for fast debugging with xargs
test_fn()
{
	echo -n "-> "; for a in "$0"; do echo -n "\"$a\" "; done; echo
	# sleep 1 # show xargs parallel mode 
}
export -f test_fn

test_fn2()
{
	echo $0; echo $1; echo $2; echo $3;
}
export -f test_fn2

# echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "test_fn"