#!/bin/bash

# compile c++ fit routine
g++ fit_charge.C `root-config --libs --cflags` -o fit_charge

#################
### CALIB FIT ###
#################

do_calib()
{	
	channel=$0 # currently not needed

	# CALIB
	f=cosmics_7to16_32to50_1325CS_100ns_isCalib0_newBL
	# f=cosmics_1to6_extra_015C_40ns_isCalib0

	./fit_charge $f
	
}
export -f do_calib

# multi process over number of channels
# (currently not needed)
# channel_list=(1 2 3 4 5 6 7 8)
channel_list=(1)

printf "%s\n" ${channel_list[@]} | xargs -n 1 -P 1 bash -c "do_calib"

echo