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

	"\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --inputFolder 13Nov2023_10kPoints_1 --runNr 1"
	"\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --inputFolder 15Nov2023_day_10kpoints_2 --runNr 2"
	"\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --inputFolder 15Nov2023_night_10kpoints_3 --runNr 3"
	"\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --inputFolder 16Nov2023_night_10kpoints_4 --runNr 4"
	"\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --inputFolder 17Nov2023_night_10kpoints_5 --runNr 5"
	"\n--nanoSec $setWindow --isCalib $setCalib --events $setEvents --inputFolder 19Nov2023_day_10kpoints_6 --runNr 6"
	
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
