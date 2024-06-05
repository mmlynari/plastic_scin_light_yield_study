import glob, os, sys
import subprocess
import csv
import argparse
import numpy as np
import ROOT as r

from multiprocessing import Pool
import pandas as pd


#######################
###### GLOBALS ########
#######################

r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(0);
r.gStyle.SetLineScalePS(2);


##################################
###### FUNCTIONS: ANALYSIS  ######
##################################


# __ Baseline Fit____________________________
# pol0 fit in range
# return bl value and reduced chi2 as list
def bl_fit(h_Wave,t1,t2):
    f_const = r.TF1("f_const","pol0",t1,t2)
    h_Wave.Fit("f_const","RNQ")

    bl_value = f_const.GetParameter(0)
    bl_red_chi2 = f_const.GetChisquare()/f_const.GetNDF()
    #print('bl_value', bl_value)

    bl_info = [bl_value,bl_red_chi2]
    return bl_info

# __ Get Maximum in Range _________________________________
# Returns amplitude value of maximum in given range t1-t2
# using a constant fit over a 0.5 ns range around the maximum.
# baseline corrected
def max_inRange(h_Wave,t1,t2,bl):
    
    f_const = r.TF1("f_const","pol0",t1,t1)

    h_Wave.GetXaxis().SetRangeUser(t1,t2)

    r1 = h_Wave.GetXaxis().GetBinCenter( h_Wave.GetMaximumBin() ) - 0.5
    r2 = r1+1

    h_Wave.Fit("f_const","RNQ","",r1,r2) # fit constant
    max = f_const.GetParameter(0) - bl # substract baseline

    ## mmlynari needs to comment back 
    #h_Wave.GetXaxis().UnZoom()

    return(max)

# __ Get Time at Maximum in Range _________________________________
# Returns time value of maximum amplitude in given range t1-t2
def t_max_inRange(h_Wave,t1,t2):

    h_Wave.GetXaxis().SetRangeUser(t1,t2)
    t_max = h_Wave.GetXaxis().GetBinCenter( h_Wave.GetMaximumBin() )
    #h_Wave.GetXaxis().UnZoom()

    return(t_max)

# __ Get Integral in Range _________________________________
# assumes baseline corrected histogram
# possible additional correction: offset
def integral_inRange(h_Wave,t1,t2,offset):

    b1 = h_Wave.FindBin(t1) 
    b2 = h_Wave.FindBin(t2)
    

    
    intgrl = h_Wave.Integral(b1,b2,"width")- offset*(t2-t1)
    # alternative:
    # integral, substract baseline, substract lower and upper half of boundary bins
    # c1 = h_Wave.GetBinContent(b1)
    # c2 = h_Wave.GetBinContent(b2)
    # intgrl = h_Wave.Integral(b1,b2,"width")- offset*(t2-t1) - c1*(t1 - h_Wave.GetXaxis().GetBinLowEdge(b1)) - c2*(h_Wave.GetXaxis().GetBinLowEdge(b2)-t2 )


    return(intgrl)

# __ Get Integral in Range _________________________________
# integral baseline correction:
# takes baseline integral and subtracts from signal integral
def integral_inRange_alt(h_Wave,t1,t2,bl_offset,bl_window):

    b1 = h_Wave.FindBin(t1) 
    b2 = h_Wave.FindBin(t2)

    # baseline integral range
    b1_alt = h_Wave.FindBin(bl_offset) 
    b2_alt = h_Wave.FindBin(bl_offset+bl_window)
    # to scale to signal integration window size
    window_ratio = (t2-t1) / bl_window

    intgrl = h_Wave.Integral(b1,b2,"width") - h_Wave.Integral(b1_alt,b2_alt,"width")*window_ratio


    return(intgrl)

# __ CONSTANT FRACTION _________________________________
# Get Signal Time
def cfd(h_Wave,thr):

    peak = h_Wave.GetMaximum()
    
    timePos = 1
    val = 0
    # time [bins] when threshold is reached
    while abs(val)<thr*peak: 
        timePos+=1
        val = h_Wave.GetBinContent(timePos)
        # print(val)

    # interpolate in between bins:
    # x1 = h_Wave.GetBinCenter(timePos-10)
    # x2 = h_Wave.GetBinCenter(timePos)
    # y1 = h_Wave.GetBinContent(timePos-10)
    # y2 = h_Wave.GetBinContent(timePos)

    # k = (x2-x1)/(y2-y1)
    # signal_time = x1+k*(thr*peak-y1)

    signal_time = h_Wave.GetBinCenter(timePos)

    return signal_time

def cfd_t_fall(h_Wave,thr):
    peak = h_Wave.GetMaximum()
    # print(peak)
    timePos = h_Wave.GetMaximumBin()
    val = peak
    while abs(val)>thr*peak:
        timePos += 20
        val = h_Wave.GetBinContent(timePos)
        # print("val = %f"%val)
        # print("time = %f"%h_Wave.GetBinCenter(timePos))
    
    signal_fall_time = h_Wave.GetBinCenter(timePos)
    return signal_fall_time


#####################################
###### FUNCTIONS: CSV TO HISTO ######
#####################################

# read csv file
# write to pandas dataframe
def csv2pd(csv_file):
    data = pd.read_csv(csv_file,delimiter=",",names = ["time", "amp_sipm", "amp_trigger"],skiprows=21)
    data["time"] = data["time"].multiply(10**9) # convert time to nanoseconds
    if len( data["amp_trigger"].dropna() ) == 0 : # trigger included in data set ?
        data = data.drop(["amp_trigger"],axis=1)
    else:
        data["amp_trigger"] = data["amp_trigger"].multiply(-1) # invert waveform
        data["amp_trigger"] = data["amp_trigger"].add(3.5) # shift to positive values
    return data

# fill histograms for each channel from dataframe
# return list of histograms
def pd2hist_list(df):
    hist_list = []
    nCh = df.shape[1]-1 # number of channels = number of columns - 1
    for j in range(0,nCh):
        # initialize histogram
        # first column for binning and range
        # second column and above filled to histogram
        h_Wave = r.TH1D("h_Wave_"+str(j+1),"h_Wave_"+str(j+1),len(df),df.iloc[0,0],df.iloc[-1,0])
        counter = 0
        for i in df.iloc[:,j+1].values:
            h_Wave.SetBinContent(counter,i)
            counter += 1
        hist_list.append(h_Wave)
    return hist_list


###########################
###### INTITIALIZE ########
###########################

directory = os.path.join( os.getcwd() )
out_directory = directory

parser = argparse.ArgumentParser()
parser.add_argument('--inputFolder', type=str, help='Specify name of the input folder with csv files.')
parser.add_argument('--events', type=int, default = 10000, help='Specify the number of events.')
parser.add_argument('--nanoSec', type=float, default=40, help='Specify the integration time.')
parser.add_argument('--isCalib', type=int, default=0, help='Specify if (1) or if not (0) charge/amplitude should be converted to N_pe.')
parser.add_argument('--runNr', type=int, default=0, help='Set the run number')


args, _ = parser.parse_known_args()

inputDir = "/eos/user/t/tilepmt/SiPM/data/csv/data_10kpoints/"+args.inputFolder 
#DC" #topTile_10mV"#noSignal" #600muV" #800muV/SiPM" #50Ohm"
print(inputDir)

root_dir = str(out_directory).replace('csv', 'root')
if not os.path.isdir(root_dir):
    os.system("mkdir -p %s"%root_dir)
outFileCombined = str(root_dir+"/ntuple_integral_"+str(int(args.nanoSec))+"ns_"+str(args.events)+"ev"+"_calib"+str(args.isCalib)+"_runNr"+str(args.runNr))+".root"
print(outFileCombined)

# directory to store waveforms
hist_dir = str(root_dir+"/waveforms")
if not os.path.isdir(hist_dir):
    os.system("mkdir -p %s"%hist_dir)
# print rates

# directory to store event trees 
event_tree_dir = str(root_dir+"/event_trees")
if not os.path.isdir(event_tree_dir):
    os.system("mkdir -p %s"%event_tree_dir)

nanoSec = args.nanoSec
print("analyze %d events"%args.events)
print('integrate over %s ns'%nanoSec)
if args.isCalib:
    print("calibrate to N_pe? YES")
else:
    print("calibrate to N_pe? NO")

# ___ SETTINGS ___
wave_print_rate = 1
print_sum = 0

is_newSiPM = 1 # 1 -> 1325CS SiPM
# is_newSiPM = 0 # 0 -> 015C SiPM

## mmlynari check these values
# calibration & pedestal correction hard coded
calib_factor_1325CS = 0.073767 # V x ns # 100 ns int. window
baseline_offset_13325 = 0.00633756 # V x ns
calib_factor_015C = 0.015343 # 20 ns int. window
baseline_offset_015C = 0.00397152 # 20 ns

if args.isCalib:
    if is_newSiPM:
        calib_factor = calib_factor_1325CS
        baseline_offset = baseline_offset_13325
        print("calibrate 1325CS")
    else:
        calib_factor = calib_factor_015C
        baseline_offset = baseline_offset_015C
        # baseline_offset = 0
        print("calibrate 015C")
else:
    calib_factor = 1
    baseline_offset = 0

# get x-axis range
try:
    temp_hist_df = csv2pd(inputDir+"/tek0003CH2.csv")
    strip_suffix = "CH2.root"
    pass
#except FileNotFoundError:
except:
    temp_hist_df = csv2pd(inputDir+"/tek0003ALL.csv")
    strip_suffix = "ALL.root"
    pass

waves_x_min = temp_hist_df.iloc[0,0]
waves_x_max = temp_hist_df.iloc[-1,0]
print("event histogram x-range: %f,%f"%(waves_x_min,waves_x_max))
del temp_hist_df


# __ cosmics run 7-16,32-39 __
amp_range_l = -500
amp_range_u = 0
dc_point = waves_x_min+5+0.25*nanoSec
bl_window = 300
## mmlynari array of the values in the range between dc_point+nanoSec and waves_x_max-bl_window spaced by bl_window
search_bl_array = np.arange(dc_point+nanoSec,waves_x_max-bl_window,bl_window)
## search_bl_array = np.arange(waves_x_max-3*bl_window,waves_x_max-2*bl_window,bl_window)

print_wf_range_l = waves_x_min
print_wf_range_u = waves_x_max

#print(print_wf_range_l,print_wf_range_u)

# trigger threshold
# for time over threshold and trigger time resolution
trig_thr = 0.3

#############################
###### ANALYZE HSITO ########
#############################

# 1. create a root file per event, multithreaded
# key branches: 
#   - SiPM signal integral (charge_alt)
#   - dark count integral (dc_charge_alt)
#   - baseline value for background reduction (bl_value)
# 2. combine all event trees to single tree

def analyze_hist(filename):

    # event tree
    tree_name = filename.replace("csv","root").split("/")[-1]
    outFilePerEvt = r.TFile( event_tree_dir+"/"+tree_name, 'recreate')
    tree = r.TTree('ntuple', 'ntuple')

    # extrct event number form tree_name
    # tree_name format: "tekddddALL.root"
    event_nr = 0
    event_nr_str = tree_name.replace("tek","").replace(strip_suffix,"")
    # print("event nr: {0}".format(event_nr))
    if event_nr_str == "0000":
        event_nr_str = event_nr_str.replace(event_nr_str[:3],"")
        event_nr = int(event_nr_str)
        print("event nr: {0}".format(event_nr))
    else:
        try:
            event_nr_str = event_nr_str.lstrip("0")
            event_nr = int(event_nr_str)
        except ValueError:
            print("ValueError at event {0}".format(event_nr))
            pass

    # initialize variables to store in branch
    run_nr = np.zeros(1, dtype=int)
    bl_value = np.zeros(1, dtype=float)
    bl_int = np.zeros(1, dtype=float)
    bl_rchi2 = np.zeros(1, dtype=float)
    charge =  np.zeros(1, dtype=float)
    charge_alt =  np.zeros(1, dtype=float)
    dc_charge =  np.zeros(1, dtype=float)
    dc_charge_alt =  np.zeros(1, dtype=float)
    max_amp =  np.zeros(1, dtype=float)
    t_max_amp = np.zeros(1, dtype=float)
    t_trig =  np.zeros(1, dtype=float)
    t_trig_fall =  np.zeros(1, dtype=float)
    trig_length =  np.zeros(1, dtype=float)

    # initialize branches
    tree.Branch('run_nr', run_nr, "run_nr/I")
    tree.Branch('bl_value', bl_value, "bl_value/D")
    tree.Branch('bl_int', bl_int, "bl_int/D")
    tree.Branch("bl_rchi2", bl_rchi2, "bl_rchi2/D")
    tree.Branch('charge_alt', charge_alt, "charge_alt/D")
    tree.Branch('charge', charge, "charge/D")
    tree.Branch('dc_charge', dc_charge, "dc_charge/D")
    tree.Branch('dc_charge_alt', dc_charge_alt, "dc_charge_alt/D")
    tree.Branch('max_amp', max_amp, "max_amp/D")
    tree.Branch('t_max_amp', t_max_amp, "t_max_amp/D")
    tree.Branch('t_trig', t_trig, "t_trig/D")
    tree.Branch('t_trig_fall', t_trig_fall, "t_trig_fall/D")
    tree.Branch('trig_length', trig_length, "trig_length/D")

    # mmlynari get run number (taken as the last number after _ in the input folder name)
    ## run_nr[0] = 1
    run_nr[0] = int(args.inputFolder.split("_")[-1])
    ## print("run_nr: ", run_nr)

    # read csv to pandas dataframe to root histogram
    df = csv2pd(filename)
    h_Wave_list = pd2hist_list(df)

    # h_Wave_list[0].GetXaxis().SetRangeUser(waves_x_min,waves_x_max)
    # t_low = h_Wave_list[0].GetBinCenter(1)
    # t_hig = h_Wave_list[0].GetBinCenter(10000)
    # print("hist range: %f - %f | event_nr : %d"%(t_low,t_hig,event_nr))

    #___ COMPUTE VARIABLES ____

    # time of signal maximum
    cloneHist = h_Wave_list[0].Clone() # clone, avoid potential change of histogram ranges
    t_max_amp[0] = t_max_inRange(cloneHist,amp_range_l,amp_range_u)
    # t_max_amp = t_max_inRange(h_Wave_list[0],amp_range_l,amp_range_u)

    # charge integration window
    # everything is relative to this window
    int_range_l = t_max_amp[0] - int(0.25*args.nanoSec)
    int_range_u = t_max_amp[0] + int(0.75*args.nanoSec) 

    # get clean baseline window, scan over waveform
    min_chi2 = 10
    offset = 0
    for i in search_bl_array:
        #print(i)
        t1 = i
        t2 = i + bl_window
        #print("event: %d | WINDOW: %1.1f %1.1f"%(event_nr,t1,t2))
        temp_list = bl_fit(h_Wave_list[0],t1,t2)        
        #print("bl = %1.4f , rchi2 = %1.6f"%(temp_list[0],temp_list[1]))

        temp_rchi2 = temp_list[1]
        if temp_rchi2 < min_chi2:
            min_chi2 = temp_rchi2
            offset = +i
            # print(offset)

    # baseline
    bl_range_l = offset
    bl_range_u = offset+bl_window
    bl_info = bl_fit(h_Wave_list[0],bl_range_l,bl_range_u)
    bl_value[0] = bl_info[0]
    bl_rchi2[0] = bl_info[1]

    # baseline correction for entire hsitogram
    # (not needed in the moment, integral baseline preferred -> integral_inRange_alt)
    # bl_shift = r.TF1("bl_shift","pol0",waves_x_min,waves_x_max)
    # bl_shift.SetParameter(0,bl_value)
    # h_Wave_list[0].Add(bl_shift,-1)

    # integral over baseline range
    bl_int[0] = ( integral_inRange(h_Wave_list[0],bl_range_l,bl_range_u,0) - baseline_offset ) / calib_factor

    # amplitude
    max_amp[0] = max_inRange(cloneHist,amp_range_l,amp_range_u,bl_value[0])
    # max_amp[0] = max_inRange(h_Wave_list[0],amp_range_l,amp_range_u,bl_value[0])

    # charge
    charge[0] = ( integral_inRange(h_Wave_list[0],int_range_l,int_range_u,0) - baseline_offset ) / calib_factor
    charge_alt[0] = ( integral_inRange_alt(h_Wave_list[0],int_range_l,int_range_u,bl_range_l,bl_window) - baseline_offset ) / calib_factor

    # dark counts
    dc_range_l = dc_point - int(0.25*args.nanoSec)
    dc_range_u = dc_point + int(0.75*args.nanoSec)
    dc_charge[0] = ( integral_inRange(h_Wave_list[0],dc_range_l,dc_range_u,0) - 0 ) / calib_factor
    dc_charge_alt[0] = ( integral_inRange_alt(h_Wave_list[0],dc_range_l,dc_range_u,bl_range_l,bl_window) - 0 ) / calib_factor
    #print("dc_point", dc_point) #360

    # trigger time
    if len(h_Wave_list)>1: 
        t_trig[0] = cfd(h_Wave_list[1],trig_thr)
        t_trig_fall[0] = cfd_t_fall(h_Wave_list[1],1-trig_thr)
        trig_length[0] = t_trig_fall[0]-t_trig[0]
        # print("t_trig = %f"%t_trig)

    #___ PRINT WAVEFORMS ____

    if event_nr == event_nr and charge_alt[0]>-999.: 
        #print("event_nr: ", event_nr)
        #print("wave_print_rate", wave_print_rate)
    # additional filter possible: 
    # if event_nr%wave_print_rate == 0 or dc_charge<-0.7:
        c_waves = r.TCanvas("c_waves","c_waves",600,500)
        c_waves.SetLeftMargin(1.4);
        c_waves.SetBottomMargin(1.4);
        c_waves.Divide(1,len(h_Wave_list))
        
        # SiPM waveform
        c_waves.cd(1)

        h_Wave_list[0].GetXaxis().SetRangeUser(print_wf_range_l,print_wf_range_u)
        h_Wave_list[0].SetTitle("cosmics measurement - event waveform - SiPM: 1325CS")
        h_Wave_list[0].GetXaxis().SetTitle("time [ns]")
        h_Wave_list[0].GetYaxis().SetTitle("amplitude[V]")

        h_Wave_list[0].DrawCopy()

        # baseline
        ln_bl = r.TLine(bl_range_l,bl_value[0],bl_range_u,bl_value[0])
        # maximum amplitude
        ln_max_amp = r.TLine(t_max_amp[0],bl_value[0],t_max_amp[0],max_amp[0])
        # charge integration window
        ln_intW_l = r.TLine(int_range_l,bl_value[0],int_range_l,max_amp[0])
        ln_intW_u = r.TLine(int_range_u,bl_value[0],int_range_u,max_amp[0])
        # dc integration window
        ln_dcW_l = r.TLine(dc_range_l,bl_value[0],dc_range_l,max_amp[0])
        ln_dcW_u = r.TLine(dc_range_u,bl_value[0],dc_range_u,max_amp[0])
        # baseline offset integration window
        ln_blOff_l = r.TLine(offset,bl_value[0],offset,max_amp[0])
        ln_blOff_u = r.TLine(bl_window+offset,bl_value[0],bl_window+offset,max_amp[0])

        ln_bl.SetLineColor(2)
        ln_bl.SetLineStyle(2)
        ln_max_amp.SetLineColor(3)
        ln_max_amp.SetLineStyle(9 )
        ln_intW_u.SetLineColor(8)
        ln_intW_l.SetLineColor(8)
        ln_dcW_u.SetLineColor(7)
        ln_dcW_l.SetLineColor(7)
        ln_dcW_u.SetLineStyle(5)
        ln_dcW_l.SetLineStyle(5)
        ln_blOff_l.SetLineColor(2)
        ln_blOff_u.SetLineColor(2)
        ln_blOff_l.SetLineStyle(2)
        ln_blOff_u.SetLineStyle(2)
        ln_list = [ln_bl,ln_max_amp,ln_intW_u,ln_intW_l,ln_dcW_u,ln_dcW_l,ln_blOff_l,ln_blOff_u]
        for line in ln_list:
            line.Draw("same")

        h_leg = r.TLegend(0.6,0.7,0.9,0.9);
        h_leg.SetTextSize(0.02);
        h_leg.AddEntry(h_Wave_list[0],r.Form("waveform data"),"l");
        h_leg.AddEntry(ln_bl,r.Form("baseline window"),"l");
        h_leg.AddEntry(ln_max_amp,r.Form("max. amplitude"),"l");
        h_leg.AddEntry(ln_intW_l,r.Form("signal window"),"l");
        h_leg.AddEntry(ln_dcW_u,r.Form("dark count window"),"l");
        h_leg.Draw()


        # trigger wafeform
        if len(h_Wave_list)>1:

            c_waves.cd(2)
            h_Wave_list[1].SetTitle("cosmics measurement - event waveform - trigger: Cosmic Pi")
            h_Wave_list[1].GetXaxis().SetTitle("time [ns]")
            h_Wave_list[1].GetYaxis().SetTitle("amplitude[V]")
            h_Wave_list[1].GetXaxis().SetRangeUser(print_wf_range_l,print_wf_range_u)
            h_Wave_list[1].DrawCopy()

            # trigger time
            ln_t_trig = r.TLine(t_trig[0],0,t_trig[0],4.5)
            ln_t_trig_fall = r.TLine(t_trig_fall[0],0,t_trig_fall[0],4.5)

            ln_t_trig.SetLineColor(3)
            ln_t_trig_fall.SetLineColor(3)

            ln_list = [ln_t_trig,ln_t_trig_fall]
            for line in ln_list:
                line.Draw("same")
        #if max_amp[0]>0.04:
        #c_waves.Print(hist_dir+"/wave_ev"+str(event_nr)+".png")

    # store in ROOT tree
    tree.Fill()

    # write event tree
    tree.Write()        
    outFilePerEvt.Write()
    outFilePerEvt.Close()

# get list of input filenames
file_list = glob.glob(inputDir+'/*.csv')

# number of parallel threads, depends on CPU cores
pool = Pool(processes=10) 
pool.map(analyze_hist,file_list[0:args.events])

######################
###### EXPORT ########
######################

# store all event trees in a combined root file
sum_hist_SiPM = r.TH1D("sum_hist_SiPM","sum_hist_SiPM",1000,waves_x_min,waves_x_max)
# sum_hist_SiPM = r.THStack("sum_hist_SiPM","sum of SiPM waveforms;time [ns]")

chain = r.TChain("ntuple")
root_file_list = glob.glob(event_tree_dir+"/*"+strip_suffix)

for i in range(0,args.events):
    chain.Add(root_file_list[i])

    # waveform sum plot not working in the moment...
    
    # if (print_sum & i<100):
    #     event_file = r.TFile.Open(root_file_list[i])
    #     temp_h = r.TH1D("temp_h","temp_h",10000,waves_x_min,waves_x_max)
    #     event_file.GetObject("h_Wave_1",temp_h)
    #     for i in range(0,10000):
    #         # new_bin_content = sum_hist_SiPM.GetBinContent(i)+temp_h.GetBinContent(i)
    #         new_bin_content =temp_h.GetBinContent(i)
    #         sum_hist_SiPM.AddBinContent(i,new_bin_content)
    #     # temp_h.GetXaxis().SetRangeUser(waves_x_min,waves_x_max)
    #     # sum_hist_SiPM.Add(temp_h) 
    #     event_file.Close()
    #     del temp_h

chain.Merge(outFileCombined)
# new_tree = r.TFile.Open(outFileCombined,"update")
# sum_hist_SiPM.Draw()
# sum_hist_SiPM.Write()
# new_tree.Close()



# if print_sum:
#     sum_C = r.TCanvas("sum_C","sum_C")
#     sum_C.cd()
#     # sum_hist_SiPM.GetStack().Last().Draw()
#     # for i in range(0,10000):
#     #     test.SetBinContent(i,sum_hist_SiPM.GetBinContent(i))
#     sum_hist_SiPM.Draw()
#     r.gPad.Update()
#     r.gPad.SaveAs(hist_dir+"/sum_wave.pdf","pdf")


# merge_command = "hadd -f -k {0} {1}/*ALL.root".format(outFileCombined,root_dir)
# # subprocess.run(merge_command)
# os.system(merge_command)
# sys.exit("done")

