########################################################################
######################### LDHELMET - PADE STEP #########################
########################################################################


#!/bin/bash

## Input file : output.conf file
## Output file : output.pade
## Define theta
tetha=0.01
## Default parameters are used of the number of pade coefficient -x 

time ./ldhelmet pade --num_threads 20 -t ${theta} -x 11 --defect_threshold 40 -c output.conf -o output.pade
