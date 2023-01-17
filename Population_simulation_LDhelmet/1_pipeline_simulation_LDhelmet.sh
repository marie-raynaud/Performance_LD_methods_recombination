########################################################################

######## PIPELINE POPULATION SIMULATION AND LDHELMET INFERENCE #########
########               CONSTANT POPULATION SIZE                #########

########################################################################

#!/bin/bash

########## DESCRIPTION ##########
### This script launch the scripts of the LDhelmet pipeline.
## It launch the different command lines and scripts to simulate populations from know underlying recombination landscapes using the coalescent program MSPRIME under various methodological and demographic parameters.
## It then launch the scripts to run the 5 LDhelmet steps to infer the recombination rates of the simulatedp populations.
## Finally, it smooth the recombination landscapes in 500 and 2500 bp windows. This step is optional but has been done to analyse the landscapes.
## The parameters Ne, SS, BP, theta, underlying_landscapes, need to be defined at the beginning of each scripts of the pipeline. 

## The data can be organized in 4 folders : 
# - scripts : containing all the scripts of the pipeline
# - input_LDhelmet : containing all the input file for LDhelmet that are generating during the pipeline
# - output_LDhelmet : containing the output files of LDhelmet
# - Underlying_landscapes : the underlying landscapes to be used of the simulation of the populations in MSPRIME

## The script can be launched in a loop to generates replicated of the simulated populations with the same set of parameters. 




########## PIPELINE ##########
### Run the simulation with Msprime and format the input files for LDhelmet by launching the command line of the script:
2_input_LDhelmet.sh 

### Create the ancestral prior states matrix for LDhelmet
Rscript 5_anc_prior.R

### Launch the scripts to run the different steps of LDhelmet, to infer the recombination landscapes of the simulated populations
## First the user needs to check, and to change, the variable used (i.e. input/output files, theta, block penalty)
bash 6_run_all_ldhelmet.bash

### Smooth recombination rates in 500 and 2500 bp windows
## This step is optional or can be modified by the user. 
# Retrieve the columns containing the positions of the SNPs (the left SNPs) and the mean recombination rates infered by LDhelmet in a new file 
grep -v '#' output.txt | grep -v 'ver' | cut -d ' ' -f 1,3 > LDhelmet_raw_rates.txt
# Run a python script to convert the underlying and simulated rates stored in the file Underlying_Simulated_Map.txt in units of rho/bp, and smooth the underlying, simulated and inferred recombination rates in windows of 500 bp and 2.5kb
python 7_convert_smooth_map.py





