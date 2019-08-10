#!/bin/bash

setCal=1
setWindow=100
setDir="/Users/julianschliwinski/Studium/CERN_Summer_2019/analysis/cosmics_piTrigger/root"
setFile="ntuple_integral_"$setWindow"ns_10000ev_calib$setCal.root"

setOutDir="/Users/julianschliwinski/Studium/CERN_Summer_2019/analysis/cosmics_piTrigger/comb_root"

##############
### 1325CS ###
##############
# __ ALL ___
outFile="$setOutDir/cosmics_7to16_32to50_1325CS_"$setWindow"ns_isCalib"$setCal"_newBL.root"
hadd $outFile $setDir/cosmics_middle_7/$setFile $setDir/cosmics_middle_8/$setFile $setDir/cosmics_middle_9/$setFile $setDir/cosmics_middle_10/$setFile $setDir/cosmics_middle_11/ntuple_integral_"$setWindow"ns_5336ev_calib$setCal.root $setDir/cosmics_middle_12/$setFile $setDir/cosmics_middle_13/$setFile $setDir/cosmics_middle_14/$setFile $setDir/cosmics_middle_15/ntuple_integral_"$setWindow"ns_9972ev_calib$setCal.root $setDir/cosmics_middle_16/ntuple_integral_"$setWindow"ns_9932ev_calib$setCal.root $setDir/cosmics_middle_32/ntuple_integral_"$setWindow"ns_7073ev_calib$setCal.root $setDir/cosmics_middle_33/ntuple_integral_"$setWindow"ns_6429ev_calib$setCal.root $setDir/cosmics_middle_34/ntuple_integral_"$setWindow"ns_3791ev_calib$setCal.root $setDir/cosmics_middle_35/ntuple_integral_"$setWindow"ns_8833ev_calib$setCal.root $setDir/cosmics_middle_36/ntuple_integral_"$setWindow"ns_5235ev_calib$setCal.root $setDir/cosmics_middle_37/ntuple_integral_"$setWindow"ns_8485ev_calib$setCal.root $setDir/cosmics_middle_38/ntuple_integral_"$setWindow"ns_7386ev_calib$setCal.root $setDir/cosmics_middle_39/ntuple_integral_"$setWindow"ns_9206ev_calib$setCal.root $setDir/cosmics_middle_40/ntuple_integral_"$setWindow"ns_5115ev_calib$setCal.root $setDir/cosmics_middle_41/ntuple_integral_"$setWindow"ns_6839ev_calib$setCal.root $setDir/cosmics_middle_42/ntuple_integral_"$setWindow"ns_8286ev_calib$setCal.root $setDir/cosmics_middle_43/ntuple_integral_"$setWindow"ns_3808ev_calib$setCal.root $setDir/cosmics_middle_44/ntuple_integral_"$setWindow"ns_6821ev_calib$setCal.root $setDir/cosmics_middle_45/ntuple_integral_"$setWindow"ns_4917ev_calib$setCal.root $setDir/cosmics_middle_46/ntuple_integral_"$setWindow"ns_9825ev_calib$setCal.root $setDir/cosmics_middle_47/ntuple_integral_"$setWindow"ns_9791ev_calib$setCal.root $setDir/cosmics_middle_48/ntuple_integral_"$setWindow"ns_5716ev_calib$setCal.root $setDir/cosmics_middle_49/ntuple_integral_"$setWindow"ns_10000ev_calib$setCal.root $setDir/cosmics_middle_50/ntuple_integral_"$setWindow"ns_9242ev_calib$setCal.root

# outFile="$setOutDir/cosmics_7to16_32to50_1325CS_"$setWindow"ns_isCalib"$setCal"_newBL.root"
# hadd $outFile $setOutDir/cosmics_7to16_32to46_1325CS_100ns_isCalib0_newBL.root $setDir/cosmics_middle_47/ntuple_integral_"$setWindow"ns_9791ev_calib$setCal.root $setDir/cosmics_middle_48/ntuple_integral_"$setWindow"ns_5716ev_calib$setCal.root $setDir/cosmics_middle_49/ntuple_integral_"$setWindow"ns_10000ev_calib$setCal.root $setDir/cosmics_middle_50/ntuple_integral_"$setWindow"ns_9242ev_calib$setCal.root


# # ___ FCC_long_10, short fiber ___
# outFile="$setOutDir/cosmics_7to8_1325CS_80ns_isCalib$setCal.root"
# hadd $outFile $setDir/cosmics_middle_7/$setFile $setDir/cosmics_middle_8/$setFile

# # ___ FCC_1_a+FCC_1_b, short fiber ___
# outFile="$setOutDir/cosmics_10to11_1325CS_80ns_isCalib$setCal.root"
# hadd $outFile $setDir/cosmics_middle_10/$setFile $setDir/cosmics_middle_11/ntuple_integral_80ns_5336ev_calib$setCal.root

# ___ FCC_1_a+FCC_1_b, long fiber ___
# outFile="$setOutDir/cosmics_12to15_1325CS_80ns_isCalib$setCal.root"
# hadd $outFile $setDir/cosmics_middle_12/$setFile $setDir/cosmics_middle_13/$setFile $setDir/cosmics_middle_14/$setFile $setDir/cosmics_middle_15/ntuple_integral_80ns_9972ev_calib$setCal.root

# # ___ FCC_long_10, Mylar, long fiber ___
# outFile="$setOutDir/cosmics_34to39_1325CS_"$setWindow"ns_isCalib"$setCal"_newBL.root"
# hadd $outFile $setDir/cosmics_middle_34/ntuple_integral_"$setWindow"ns_3791ev_calib$setCal.root $setDir/cosmics_middle_35/ntuple_integral_"$setWindow"ns_8833ev_calib$setCal.root $setDir/cosmics_middle_36/ntuple_integral_"$setWindow"ns_5235ev_calib$setCal.root $setDir/cosmics_middle_37/ntuple_integral_"$setWindow"ns_8485ev_calib$setCal.root $setDir/cosmics_middle_38/ntuple_integral_"$setWindow"ns_7386ev_calib$setCal.root $setDir/cosmics_middle_39/ntuple_integral_"$setWindow"ns_9206ev_calib$setCal.root

# ___ FCC_long_10, Alu, short fiber ___
# outFile="$setOutDir/cosmics_40to41_1325CS_"$setWindow"ns_isCalib$setCal.root"
# hadd $outFile $setDir/cosmics_middle_40/ntuple_integral_"$setWindow"ns_5115ev_calib$setCal.root $setDir/cosmics_middle_41/ntuple_integral_"$setWindow"ns_6839ev_calib$setCal.root

# ___ setup, 1c) FCC_1, Alu, long fiber ___
# outFile="$setOutDir/cosmics_47to50_1325CS_"$setWindow"ns_isCalib"$setCal"_newBL.root"
# hadd $outFile $setDir/cosmics_middle_47/ntuple_integral_"$setWindow"ns_9791ev_calib$setCal.root $setDir/cosmics_middle_48/ntuple_integral_"$setWindow"ns_5716ev_calib$setCal.root $setDir/cosmics_middle_49/ntuple_integral_"$setWindow"ns_10000ev_calib$setCal.root $setDir/cosmics_middle_50/ntuple_integral_"$setWindow"ns_9242ev_calib$setCal.root

############
### 015C ###
############

# ___ ALL ___
# outFile="$setOutDir/cosmics_1to6_015C_"$setWindow"ns_isCalib"$setCal"_newBL.root"
# hadd $outFile $setDir/cosmics_bottom_1/ntuple_integral_"$setWindow"ns_2001ev_calib$setCal.root $setDir/cosmics_bottom_2/ntuple_integral_"$setWindow"ns_2981ev_calib$setCal.root $setDir/cosmics_bottom_3/$setFile $setDir/cosmics_middle_4/$setFile $setDir/cosmics_middle_5/ntuple_integral_"$setWindow"ns_5545ev_calib$setCal.root $setDir/cosmics_middle_6/$setFile

# ___ EXTRA RUNS ___
# outFile="$setOutDir/cosmics_103to109_015C_"$setWindow"ns_isCalib$setCal.root"
# hadd $outFile $setDir/muons_tile10_2mFibre_103/ntuple_integral_"$setWindow"ns_9937ev_calib$setCal.root $setDir/muons_tile10_2mFibre_104/ntuple_integral_"$setWindow"ns_10000ev_calib$setCal.root $setDir/muons_tile10_2mFibre_105/ntuple_integral_"$setWindow"ns_7975ev_calib$setCal.root $setDir/muons_tile10_05mFibre_107/ntuple_integral_"$setWindow"ns_9999ev_calib$setCal.root $setDir/muons_tile10_05mFibre_108/ntuple_integral_"$setWindow"ns_9631ev_calib$setCal.root $setDir/muons_tile10_05mFibre_109/ntuple_integral_"$setWindow"ns_8520ev_calib$setCal.root

# ___ 1-6 + EXTRA RUNS ___
# outFile="$setOutDir/cosmics_1to6_extra_015C_"$setWindow"ns_isCalib$setCal.root"
# hadd $outFile $setOutDir/cosmics_1to6_015C_40ns_isCalib$setCal.root $setOutDir/cosmics_103to109_015C_"$setWindow"ns_isCalib$setCal.root
