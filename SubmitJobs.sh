#!/bin/bash

files_to_load=(
	# EG Data
        crab_DoubleEG_RunB           	# 0
        crab_DoubleEG_RunC           	# 1
        crab_DoubleEG_RunD           	# 2
        crab_DoubleEG_RunE           	# 3
        crab_DoubleEG_RunF           	# 4
        crab_DoubleEG_RunG           	# 5
        crab_DoubleEG_RunHver2       	# 6
        crab_DoubleEG_RunHver3       	# 7

        # MC Signal
        DYLL_M10to50              	# 8
        DYLL_M50toInf             	# 9
        DYLL_M100to200            	# 10
        DYLL_M200to400            	# 11
        DYLL_M400to500            	# 12
        DYLL_M500to700            	# 13
        DYLL_M700to800            	# 14
        DYLL_M800to1000           	# 15
        DYLL_M1000to1500          	# 16
        DYLL_M1500to2000          	# 17
        DYLL_M2000to3000          	# 18

        # Tops
        ST_tW                        	# 19
        ST_tbarW                     	# 20
	ttbar				# 21
        ttbar_M700to1000             	# 23
        ttbar_M1000toInf             	# 24

        # EW
        WW                           	# 25
        WZ                           	# 26
        ZZ                           	# 27

        # Fakes
        WJetsToLNu_amcatnlo          	# 39   
        WJetsToLNu_amcatnlo_ext      	# 40   
        WJetsToLNu_amcatnlo_ext2v5   	# 41   

	# EM Enriched QCD			
	QCDEMEnriched_Pt20to30		# 42
	QCDEMEnriched_Pt30to50		# 43
	QCDEMEnriched_Pt50to80		# 44
	QCDEMEnriched_Pt50to80_ext1	# 45
	QCDEMEnriched_Pt80to120		# 46
	QCDEMEnriched_Pt80to120_ext1	# 47
	QCDEMEnriched_Pt120to170	# 48
	QCDEMEnriched_Pt120to170_ext1	# 49
	QCDEMEnriched_Pt170to300	# 50
	QCDEMEnriched_Pt300toInf	# 51
)
txt_directory="fileLists/"
extension=".txt"

for index in ${!files_to_load[*]}; do
	echo "Beginning to process ${files_to_load[$index]}"
	echo "Loading text file ${txt_directory}${files_to_load[$index]}${extension}"
	condor_submit \
		arg1=${files_to_load[$index]} \
		arg2=${txt_directory}${files_to_load[$index]}${extension} \
		condor_control.condor
done
