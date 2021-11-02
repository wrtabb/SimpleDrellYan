#!/bin/bash

root -l << EOF
.L analyzeData.C+
analyzeData("DYLL_M500to700_TauTau")
EOF
