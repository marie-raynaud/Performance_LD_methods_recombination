########################################################################
###################### LDHELMET - FIND_CONFS STEP ######################
########################################################################

#!/bin/bash

## Input file : fasta file of the samples
## Output file : output.conf
## Window size -w and number of threads can be changed

time ./ldhelmet find_confs --num_threads 20 -w 50 -o output.conf Constant_size_pop.fa
