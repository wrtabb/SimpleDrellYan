#!/bin/bash

files_to_load=(
        crab_DoubleEG_RunB           # 0
        crab_DoubleEG_RunC           # 1
        crab_DoubleEG_RunD           # 2
        crab_DoubleEG_RunE           # 3
        crab_DoubleEG_RunF           # 4
        crab_DoubleEG_RunG           # 5
        crab_DoubleEG_RunHver2       # 6
        crab_DoubleEG_RunHver3       # 7

        # MC Signal
        DYLL_M10to50_EE              # 8
        DYLL_M50to100_EE             # 9
        DYLL_M100to200_EE            # 10
        DYLL_M200to400_EE            # 11
        DYLL_M400to500_EE            # 12
        DYLL_M500to700_EE            # 13
        DYLL_M700to800_EE            # 14
        DYLL_M800to1000_EE           # 15
        DYLL_M1000to1500_EE          # 16
        DYLL_M1500to2000_EE          # 17
        DYLL_M2000to3000_EE          # 18

        # Tops
        ST_tW                        # 19
        ST_tbarW                     # 20
        ttbar_M0to700      	     # 21
        ttbar_M700to1000             # 22
        ttbar_M1000toInf             # 23

        # EW
        WW                           # 24
        WZ                           # 25
        ZZ                           # 26
        DYLL_M10to50_TauTau          # 27
        DYLL_M50to100_TauTau         # 28
        DYLL_M100to200_TauTau        # 29
        DYLL_M200to400_TauTau        # 30
        DYLL_M400to500_TauTau        # 31
        DYLL_M500to700_TauTau        # 32
        DYLL_M700to800_TauTau        # 33
        DYLL_M800to1000_TauTau       # 34
        DYLL_M1000to1500_TauTau      # 35
        DYLL_M1500to2000_TauTau      # 36
        DYLL_M2000to3000_TauTau      # 37

        # Fakes
        WJetsToLNu_amcatnlo          # 38   
        WJetsToLNu_amcatnlo_ext      # 39   
)

for index in ${!files_to_load[*]}; do
	echo "Beginning to process ${files_to_load[$index]}"
	condor_submit \
		arg1=${files_to_load[$index]} \
		condor_control.condor
done
