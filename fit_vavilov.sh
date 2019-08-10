#!/bin/bash

# compile c++ fit routine
g++ fit_vavilov.C `root-config --libs --cflags` -lMathMore -o fit_vavilov

#################
### CALIB FIT ###
#################

do_vavilov_fit()
{	
	filename=$0

	# ./fit_vavilov $filename
	./fit_vavilov
	
}
export -f do_vavilov_fit

data_list=(
	
	# cosmics

	"dummy_filename"
	# "\ndummy_filename2",
	# "\ndummy_filename3"
	)

# hand over list of to be unalized data sets 
# (currently not needed, data sets are defined in .C file)
echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "do_vavilov_fit"

echo