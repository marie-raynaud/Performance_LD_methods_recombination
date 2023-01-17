########################################################################
######################## ANCESTRAL STATES PRIOR ########################
########################################################################

########## DESCRIPTION ##########
### This script generate the ancestral allele matrix needed for LDhelmet which contains the probability of the 4 alleles A, C, G, T of being ancestral at the variable positions. 
## As the data simulated are perfectly polarized, the ancestral allele recorded in the VCF is set to have a probability of 1 of being ancestral. 


########## CREATE ANCESTRAL STATES MATRIX ##########

### Import file containing the positions of the VCF and the corresponding ancestral allele (POS.txt file)
pos <- read.table("./input_LDhelmet/POS.txt", header = FALSE, sep="\t")

### Generate the ùatrix containing the probability of each allele to be ancestral

anc <- NULL

for (i in 1:nrow(pos)) {
  if (pos[i,2] == "A")
  {anc <- rbind(anc, c(pos[i,1], c("1.0", "0.0", "0.0", "0.0")))}
  if (pos[i,2] == "C")
  {anc <- rbind(anc, c(pos[i,1], c("0.0", "1.0", "0.0", "0.0")))}
  if (pos[i,2] == "G")
  {anc <- rbind(anc, c(pos[i,1], c("0.0", "0.0", "1.0", "0.0")))}
  if (pos[i,2] == "T")
  {anc <- rbind(anc, c(pos[i,1], c("0.0", "0.0", "0.0", "1.0")))}
}


### Store the matrix in the text file anc_prior.txt
write.table(anc,file="./input_LDhelmet/anc_prior.txt", sep =" ", quote = F, dec = ".", col.names = FALSE, row.names = FALSE)
