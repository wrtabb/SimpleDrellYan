#!/bin/bash

root -l << EOF
.L analyzeData.C+
analyzeData("DYLL_M1000to1500_EE")
EOF
