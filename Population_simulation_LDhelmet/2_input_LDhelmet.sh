###############################################################################
############################# MAKE LDHELMET FILE ##############################
###############################################################################

#!/bin/bash

########## DESCRIPTION ##########
### This script run the script to simulate population samples with a known recombination rate under various parameters using MSPRIME.
## From the output VCF, LDhelmet input files are prepared
## The following script is made for constant population size simulations.



### Simulate population samples using MSPRIME
# Run the python script 3a_msprime_constant_size.py which simulate a population sample with a Ne, a sample size SS, with a mutation rate mu and under a know underlying landscape
# This script simulate a population sample of effecitve size Ne = 250000, sample size SS = 20, mutation rate = 10^-8, the underling landscape 1, and with a constant population size. 
python 3a_msprime_constant_size.py
# The following step are adapted to these parameters (e.g. theta = 4*Ne*u = 0.01 in the following steps)
# This script can be modified to simulate a constant size population with different Ne, SS, recombination rate and mutation rate parameters. 
# To simulate a population undergoing a bottleneck or an admixture event, use the example scripts 3b_msprime_bottleneck.py and 3c_msprime_admixture.py.

### Remove duplicated lines in the VCF
grep '#' Constant_size_pop_tmp.vcf > Constant_size_pop.vcf # keep the header of the vcf
grep -v '#' Constant_size_pop_tmp.vcf |sort -bunk2,2 >> Constant_size_pop.vcf


### Create fasta sequences for each samples of the simulated population from the VCF for LDhelmet
# Retrieve positions and reference alleles in a POS.txt file
grep -v '#' Constant_size_pop.vcf | cut -f 2,4 > POS.txt 

# Create a reference genome 
# To do that, create a 1 Mb sequence of "A" nucleotides and replace the variant positions by the reference allele, by running the script 4_make_RefGenome.py
python 4_make_RefGenome.py

# Create a fasta sequence for each sample of the population from the VCF file and the reference gneome
# Convert the vcf file into a fasta file per sample using the vcf2fasta function from vcflib.
vcf2fasta Constant_size_pop.vcf -f REF.fa

# Concatenate the fasta of all the samples to generate one fasta file containing the fasta sequences of the SS samples
cat tsk* > Constant_size_pop.fa

# Remove indiviual fasta files
rm tsk*


# Create LDhelmet input file of the rjmcmc step
# The rjmcmc step of LDhelmet need ine input file containg the positions of the vcf, and one fasta file containing the alleles of the samples
# This two files can be generated using the --ldhelmet option of vcftools
# The command output a .pos and a .snps file
vcftools --vcf Constant_size_pop.vcf --ldhelmet --chr 1 --out Constant_size_pop




