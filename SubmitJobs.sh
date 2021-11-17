#!/bin/bash

files_to_load=(
	# Data
        SingleMuon_Run2016B		     # 0
        SingleMuon_Run2016C		     # 1
        SingleMuon_Run2016D		     # 2
        SingleMuon_Run2016E		     # 3
        SingleMuon_Run2016F		     # 4
        SingleMuon_Run2016G		     # 5
        SingleMuon_Run2016Hver2		     # 6
        SingleMuon_Run2016Hver3		     # 7		

        # MC Signal
        DYLL_M10to50_MuMu		     # 8
        DYLL_M50to100_MuMu		     # 9
        DYLL_M100to200_MuMu		     # 10
        DYLL_M200to400_MuMu		     # 11
        DYLL_M400to500_MuMu		     # 12
        DYLL_M500to700_MuMu		     # 13
        DYLL_M700to800_MuMu		     # 14
        DYLL_M800to1000_MuMu		     # 15
        DYLL_M1000to1500_MuMu		     # 16
        DYLL_M1500to2000_MuMu		     # 17
        DYLL_M2000to3000_MuMu		     # 18

        # Tops
        ST_tW				     # 19
        ST_tbarW			     # 20
        ttbar_M0to700		     	     # 21
        ttbar_M700to1000		     # 22
        ttbar_M1000toInf		     # 23

        # EW
        WW				     # 24
        WZ				     # 25
        ZZ				     # 26
        DYLL_M10to50_TauTau		     # 27
        DYLL_M50to100_TauTau		     # 28
        DYLL_M100to200_TauTau		     # 29
        DYLL_M200to400_TauTau		     # 30
        DYLL_M400to500_TauTau		     # 31
        DYLL_M500to700_TauTau		     # 32
        DYLL_M700to800_TauTau		     # 33
        DYLL_M800to1000_TauTau		     # 34
        DYLL_M1000to1500_TauTau		     # 35
        DYLL_M1500to2000_TauTau		     # 36
        DYLL_M2000to3000_TauTau		     # 37

        # Fakes
        WJetsToLNu_amcatnlo		     # 38   
	WJetsToLNu_amcatnlo_ext		     # 39
	WJetsToLNu_amcatnlo_ext2v5	     # 39

	# QCD
	QCDMuEnriched_Pt15to20		     # 41
        QCDMuEnriched_Pt20to30		     # 42
        QCDMuEnriched_Pt30to50		     # 43
        QCDMuEnriched_Pt50to80		     # 44
        QCDMuEnriched_Pt80to120	 	     # 45
        QCDMuEnriched_Pt80to120_ext1	     # 46
        QCDMuEnriched_Pt120to170	     # 47
        QCDMuEnriched_Pt120to170_backup      # 48
        QCDMuEnriched_Pt170to300             # 49
        QCDMuEnriched_Pt170to300_backup      # 50
        QCDMuEnriched_Pt170to300_ext1        # 51
        QCDMuEnriched_Pt300to470             # 52
        QCDMuEnriched_Pt300to470_ext1        # 53
        QCDMuEnriched_Pt300to470_ext2        # 54
        QCDMuEnriched_Pt470to600             # 55
        QCDMuEnriched_Pt600to800             # 56
        QCDMuEnriched_Pt600to800_backup      # 57
        QCDMuEnriched_Pt600to800_ext1        # 58
        QCDMuEnriched_Pt800to1000            # 59
        QCDMuEnriched_Pt800to1000_ext1       # 60
        QCDMuEnriched_Pt800to1000_ext2       # 61
        QCDMuEnriched_Pt1000toInf            # 62
        QCDMuEnriched_Pt1000toInf_ext1       # 63
)

for index in ${!files_to_load[*]}; do
	echo "Beginning to process ${files_to_load[$index]}"
	condor_submit \
		arg1=${files_to_load[$index]} \
		condor_control.condor
done
