###############################################################################
###############         CREATE UNDERLYING LANDSCAPE         ###################
###############                  Gamma Law                      ###############  
###############################################################################


########## DESCRIPTION ############
#### This script create the underlying landscape using a gamma law. 
## The shape and rate parameters have been defined using the ChipSeq DMC1 data from Pratto et al. (2014).
## The underlying landscape has a resolution at the 500 bp level

########## IMPORT LIBRARIES ############
library(ggplot2)

### Initiate map of 1 Mb (=2000 windows of 500 bp)
# Create a vector chromosome with the name of the chromosome of the nucleotidic segment
Chromosome <- rep("chrX", 2001)
# Create a vector of positions : a recombination rate value will be sorted in every 500 bp positions
Position_bp <- seq(0, 1000000, by = 500)  
# Create an empty vector that will contains the recombination values drawn from the gamma law
Rate <- rep(NA,2001)

# Convert the 3 vectors into a dataframe
tab <- data.frame(Chromosome, Position_bp, Rate)


### Random draw of recombination rate values along the chromosome
## To model large-scale variation of recombination rates, recombination rates of the first half of the chromosomal segment are drawn in a gamma low of mean 1 cM/Mb, and the second half in a gamma low of mean 3 cM/Mb
# Centrochromomal region with a background recombination rates of 1 cM/Mb
for (i in 1:1000){
  tab[i,3] <-rgamma(1, shape = 0.1328031, rate = 0.1326721) 
}

# Subtelomeric region with a background recombination rates of 3 cM/Mb
for (i in 1001:2000){
  tab[i,3] <-rgamma(1, shape = 0.1597632, rate = 0.05324958)
}

# The last positions must have a rate of 0
tab[2001,3]<-0
tab <- as.data.frame(tab)
# Remove lines with a rate < 0
tab <- tab[tab$Rate >= 0, ]

###  The underlying landscape can be plotted using the following code: 

plot <- ggplot(tab, aes(x=Position_bp, y=Rate))+
  geom_line(size=0.5, color="dodgerblue3") + # geom_point(size=0.5, aes(color= cumul_rate)) = permet de colorer en fonction des valeurs de la colonnes
  theme_classic() + #encadrer le graphique
  ggtitle("Underling recombination landscape")+
  ylab("Recombination Rate (M/pb)")+
  xlab("Positions (bp)")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) #pour centrer les titres des axes et du graphe
plot



### Export underlying landscape in a text file
write.table(tab,file="~/Documents/Simulation_LDhelmet/Underlying_landscapes/Landscape_1", sep ="\t", quote = F, dec = ".", col.names = TRUE, row.names = FALSE)
