########################################################################
########################### LDHELMET PIPELINE ##########################
########################################################################


#!/bin/bash


########## DESCRIPTION ##########
### This script launch all the scripts of LDhelmet (see the user manual of LDhelmet) 
## To launch the scripts, 4 input files are required : the fasta file of all the samples, the .pos and .snps files, and the file containing the matrix of the ancestral prior states. 
## For the following scripts, the following variables need to be defined : theta, block penalty
## The final output file contains the recombination rates (in unit of rho/bp) between each SNPs of the VCF file



########## LAUNCH LDHELMET PIPELINE ##########


echo "Run find_confs"
bash find_confs.bash   &&
echo "Run table_gen"
bash table_gen.bash    &&
echo "Run pade"
bash pade.bash         &&
echo "Run rjmcmc"
bash rjmcmc.bash       &&
echo "Run post_to_text"
bash post_to_text.bash 

