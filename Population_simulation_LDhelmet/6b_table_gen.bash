########################################################################
###################### LDHELMET - TABLE_GEN STEP #######################
########################################################################

#!/bin/bash

## Input file : output.conf file
## Output file : output.lk
## Define theta 
tetha=0.01
## The default parameters are used for -r (see LDhelmet user manual)

time ./ldhelmet table_gen --num_threads 20 -t ${theta} -r 0.0 0.1 10.0 1.0 100.0 -c output.conf -o output.lk 
