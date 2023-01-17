########################################################################
######################## LDHELMET - RJMCMC STEP ########################
########################################################################

#!/bin/bash

## Input file : output.lk, output.pade, .pos and .snps files, anc_prior.txt
## Output file : output.post
## Define block penalty (here is 5 to generate a fine scale landscape, use values like 100 to generate borad-scale landscape)
bp=5
## Define seed to run replicates
SEED=$(date +"%N"); echo $SEED
## Default parameters are used for the number of iterations of the rjmcmc and for the burn-in 


time ./ldhelmet rjmcmc --seed $SEED --num_threads 20 -l output.lk -p output.pade --snps_file Constant_size_pop.ldhelmet.snps --pos_file Constant_size_pop.ldhelmet.ldhelmet.pos -a anc_prior.txt -b ${bp} --burn_in 100000 -n 1000000 -o output.post
