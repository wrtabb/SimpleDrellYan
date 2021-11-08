#!/bin/bash

root -l << EOF
.L analyzeData.C+
analyzeData("DYLL_M100to200_MuMu")
EOF
