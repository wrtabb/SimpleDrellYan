#!/bin/bash

files_to_load=(
        crab_DoubleEG_RunB           
        crab_DoubleEG_RunC           
        crab_DoubleEG_RunD           
        crab_DoubleEG_RunE           
        crab_DoubleEG_RunF           
        crab_DoubleEG_RunG           
        crab_DoubleEG_RunHver2       
        crab_DoubleEG_RunHver3       
        DYLL_M10to50_EE
        DYLL_M50to100_EE
        DYLL_M100to200_EE
        DYLL_M200to400_EE
        DYLL_M400to500_EE
        DYLL_M500to700_EE
        DYLL_M700to800_EE
        DYLL_M800to1000_EE
        DYLL_M1000to1500_EE
        DYLL_M1500to2000_EE
        DYLL_M2000to3000_EE
)

for index in ${!files_to_load[*]}; do
	echo "Beginning to process ${files_to_load[$index]}"
	condor_submit \
		arg1=${files_to_load[$index]} \
		condor_control.condor
done
