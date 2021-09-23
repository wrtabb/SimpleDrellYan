#!/bin/bash

files_to_load=(
	# Data
        SingleMuon_Run2016B          # 0
        SingleMuon_Run2016C          # 1
        SingleMuon_Run2016D          # 2
        SingleMuon_Run2016E          # 3
        SingleMuon_Run2016F          # 4
        SingleMuon_Run2016G          # 5
        SingleMuon_Run2016Hver2      # 6
        SingleMuon_Run2016Hver3      # 7		

        # MC Signal
        DYLL_M10to50_MuMu              # 8
        DYLL_M50to100_MuMu             # 9
        DYLL_M100to200_MuMu            # 10
        DYLL_M200to400_MuMu            # 11
        DYLL_M400to500_MuMu            # 12
        DYLL_M500to700_MuMu            # 13
        DYLL_M700to800_MuMu            # 14
        DYLL_M800to1000_MuMu           # 15
        DYLL_M1000to1500_MuMu          # 16
        DYLL_M1500to2000_MuMu          # 17
        DYLL_M2000to3000_MuMu          # 18

        # Tops
        ST_tW                        # 19
        ST_tbarW                     # 20
        ttbar_truncated_M0To700      # 21
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
)

for index in ${!files_to_load[*]}; do
	echo "Beginning to process ${files_to_load[$index]}"
	condor_submit \
		arg1=${files_to_load[$index]} \
		condor_control.condor
done
