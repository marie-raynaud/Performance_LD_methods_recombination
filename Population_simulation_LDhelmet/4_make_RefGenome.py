########################################################################
######################## CREATE REFERENCE GENOME #######################
########################################################################

########## DESCRIPTION ##########
### This script create a reference genome from the reference allele of the VCF file of the population sample.
### It outputs a REF.fa fasta file, which contains the reference sequences with non-variable positions of the VCF being a nucleotide "A", and the variable positions being the reference allele of the VCF.



########## LOAD PACKAGES ##########
from pyfaidx import Fasta


########## CREATE REFFERENCE GENOME ##########

# Generate a reference sequence containing 1 Mb of nucleotide "A" (it could be nucleotide "C", "G," or "T") 
seq = "A" * 10**6  #a modifier si on change la taille du segment simulé
# Store the referene seuqnece in the REF.fa file in a fasta format
ref = open("REF.fa", "w")
ref.write(">" + '1' + "\n" + seq + "\n")
ref.close

# Replace the variable positions of the VCF file by the reference alleles stored in the POS.txt file 
with open('./POS.txt') as mut_table:
	with Fasta('./REF.fa', mutable=True) as fasta:
		for record in fasta:
			for line in mut_table:
				pos, base = line.rstrip().split()
				record[int(pos)-1] = base
