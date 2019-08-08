# plastic_scin_light_yield_study
light yield response of scintillating plastic tiles using measurements of cosmic muons

**recommended directory structure:**
--working_directory  
	|  
	|--csv  
	|--root  
	|--comb_root  
		|  
		|--ph_spectrum  
		|--vavilov  

**first steps:**
- raw data format: One .csv file per event
	- 1st column: time
	- 2nd column: SiPM channel
	- 3rd column: trigger channel (optional)
	- used DAQ: Tektronix MDO3025
	- store raw data in working_directory/csv
- raw data handling:
	- convert .csv to pandas dataframe
	- convert dataframe to ROOT histogram
	- analyze histogram and store computed observables in ROOT file
	- ROOT files are stored in working_directory/root
- execute all software tools from working directory 

**work flow:**
1. check for potential corrupted .csv files using get_bad_events.py
2. create run list using run_list.sh -> executes next step
2. read in all events (.csv files) and create ROOT file with observables using convert_hist_m_thread.py. this is done using python multi-threading package subprocess -> 1 thread per event
3. combine multiple runs to single ROOT file using combine_tree.sh (identify runs using run_nr variable)
4. calibrate SiPM pulse-height spectrum (signal integral -> variable: charge_alt) using fit_charge.C and fit_charge.sh
5. re-run convert_hist_m_thread.py with calibration values
6. extract most probable light yield value (MVP) from data sets using fit_vavilov.C and fit_vavilov.sh
8. misc. tools: 
	8.1 plot calibrated dark count and SiPM signal integral distributions -> show_charge.py and show_DC.py
	8.2 plot single pulse-height spectrum: plot_ph.C (execute with interactive ROOT)
	8.3 plot color coded overlay of all waveforms of a measurement run: show_all_WF.py