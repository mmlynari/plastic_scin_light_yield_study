# plastic_scin_light_yield_study
light yield response of scintillating plastic tiles using measurements of cosmic muons

**SiPM data are stored in** 
-- /eos/user/t/tilepmt/SiPM/data/csv/data_10kpoints

**recommended directory structure:**  
--working_directory 
----|  
----|--csv  
----|--root  
----|--comb_root  
--------|  
--------|--ph_spectrum  
--------|--vavilov  

**first steps:**
- raw data format: One .csv file per event
	- 1st column: time
	- 2nd column: SiPM channel
	- 3rd column: trigger channel (optional)
	- used DAQ: Tektronix MDO3025
	- store raw data in working_directory/csv
- raw data handling:
	- convert .csv to pandas dataframe
	- convert dataframe to ROOT histograms
	- analyze histograms and store computed observables in a ROOT file
	- ROOT files are stored in working_directory/root
- execute all software tools from the working directory
- NOTE: on lxplus need to install pandas locally --> pip3 install pandas //user

**workflow:**
1. get_bad_events.py --> check for potentially corrupted .csv files
2. run_list.sh --> create a run list using + executes next step
2. convert_hist_m_thread.py --> read in all events (.csv files) and create a ROOT file with observables using. this is done using python multi-threading package subprocess -> 1 thread per event
3a. NOT TO BE USED FOR THE MOMENT combine_tree.sh --> combine multiple runs to single ROOT file (identify runs using run_nr variable), store new ROOT files in working_directory/comb_root
3b. TO BE USED combine multiple output root files into a single root file manually using hadd, store new ROOT files in working_directory/comb_root
4. calibrate SiPM pulse-height spectrum (signal integral -> variable: charge_alt) using fit_charge.C and fit_charge.sh, store in working_directory/comb_root/ph_spectrum
5. re-run convert_hist_m_thread.py with calibration values, re-run combine_tree.sh
6. extract the most probable light yield value (MVP) from data sets using fit_vavilov.C and fit_vavilov.sh, store in working_directory/comb_root/vavilov
8. misc. tools:   
	8.1 plot calibrated dark count and SiPM signal integral distributions -> show_charge.py and show_DC.py  
	8.2 plot single pulse-height spectrum: plot_ph.C (execute with interactive ROOT)  
	8.3 plot color coded overlay of all waveforms of a measurement run: show_all_WF.py

**output variables in ntuple**
run_nr --> Number of run, taken as the last number after _ in the input folder
bl_value --> value of the pedestal fit
bl_int -->
bl_rchi2 --> chi squared value corresponding to the best constant fit of the baseline (signal pedestal, defined as the fit with the smallest chi squared within the sliding window)
charge_alt --> 
charge --> total charge in the muon pulse, defined as the integral over the 
dc_charge -->
dc_charge_alt -->
max_amp --> amplitude of the pulse, defined as the result of a constant fit over +- 0.5 ns around the center of the highest bin
t_max_amp --> time position of the signal maximum defined as the center of the highest bin found within a preset range
t_trig -->
t_trig_fall -->
trig_length -->
