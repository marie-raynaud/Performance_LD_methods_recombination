########################################################################
######## SIMULATION OF CONSTANT SIZE POPULATIONS WITH MSPRIME ##########
########################################################################

########## DESCRIPTION ##########
#### This script use the package MSPRIME to simulate a population sample of size SS, with a constant size Ne, using a mutation rate mu, and under a known underlying landscape, undergoing an ADMIXTURE event. 
### It first load the underlying landscape.
### It simulates a tree sequence of a population sample with an effective size Ne, a sample size SS, under the underlying landscape loaded before, underoging a admixture event for which the parameters can be defined by the user.
### It adds mutations at a rate mu. 
### This script then convert the tree sequence with the number of recombination events into a recombination map with the recombination rates in M/bp. 
### It output a VCF file containing the genotypes of the variant positions of the population sample.


########## LOAD PACKAGES ##########
import msprime
import scipy.stats
import numpy as np
import os 

########## VARIABLES ##########
## The following variables allow to simulate a population size of Ne = 250000, SS = 20 individuals, with a mutation rate mu = 10^-8, under the underlying landscape called "Landscape_1" undergoing an admixture event
## Population parameters
underlying_map="./Underlying_landscapes/Landscape_1" # path/to/underlying_landscape.txt
Ne=250000 # effective size of the population in the present
SS=20 # sample size of the population sample
mu=10**-8 # mutation rate

## Admixture parameters
## The following variables allow to simulate an admixture event that happened 500 generations ago (time_admixture), that last 1000 generations (time_admixture) between two populations that splitted 10*Ne generation ago (time_split)
time_admixture=500 #time of the admixture between pop A and B, creating the admix pop
time_split=10*Ne #time of the split of the ancestral pop into the pop A and B




########## SIMULATION ##########

### Import recombination map
rate_map = msprime.RecombinationMap.read_hapmap(map)


### Population configuration
pop_admixture = [msprime.PopulationConfiguration(sample_size=SS*2, initial_size=Ne), 
	msprime.PopulationConfiguration(sample_size=0, initial_size=Ne)] 
## At present (time = 0), 20 individuals sampled in pop 0 and 0 sampled in pop 1
# SS*2 to have the number of haplotypes

### Demographic event : admixture
## Admixture event
admixture_event = [msprime.MassMigration(time_admixture, source=0, dest=1, proportion=0.5), msprime.MassMigration(time_split, source=1, dest=0, proportion=1)]
## First MassMigration event: at time 500 (generations), admixture event where each of the lineages of pop 0 has a probability of 0.5 of moving to pop 1 = admixture event
## Second MassMigration event : at time 10*Ne, split event resulting in lineage migration of pop 1 to 0



# Simulation of the tree sequence
ts = msprime.simulate(
	Ne=Ne,
	population_configurations=pop_admixture, # include the Ne and SS parameters
	demographic_events=admixture_event,
	recombination_map=rate_map)
 

# Retrieve recombination breakpoints positions
breakpoints = np.array(list(ts.breakpoints())) #récupère les breakpoints

   
### Add mutations to the tree sequence
ts = msprime.mutate(
	ts, rate=10**-8,  # mutation rate
	model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))
	
	


### Convert the tree sequence into a recombination map
positions = np.array(rate_map.get_positions()[:-1]) # retrieve positions from the underlying landscape
rates = np.array(rate_map.get_rates()[:-1]) # retrieve recombination rates (in cM/Mb) from the underlying landscape
num_bins = 2000 # define number of bins (e.g. 2000 bins to obtains bins of 500 bp in a chromosomal segment of 1Mb)


bin_rates, bin_edges, _ = scipy.stats.binned_statistic(
    positions, rates, bins=num_bins, range=(0,1000000))

# Remove missing data from the binned positions (bin_edges) and rates (bin_rates)
x = bin_edges[:-1][np.logical_not(np.isnan(v))]
y = bin_rates[np.logical_not(np.isnan(bin_rates))]

# Convert number of recombination events (=number of breakpoints) into recombination rate (in M/bp) 
# To do that, the recombination rate can be retrieved from the number of segregating recombination events, as theta (4Neu) can be retrieved from the number of segregating sites

bin_breakpoints, bin_edges = np.histogram(breakpoints, num_bins, density=None) # get number of breakpoints in the 500 bp bins
n = 2*SS # sample_size, need to give the number of haplotypes (so 2*SS for diploid individuals)

## Extract r (or rho) from the number of breakpoints
# Getting the number of recombination breakpoints is equivalent to getting the number of segregating sites
# E(R) = breakpoints = Somme(i=n, n-1) 1/i * 4Ner * L
# Somme(i=n, n-1) 1/i = a
# L = length of the bins
# => Rho = breakpoints / a
# => r = rho/(4Ne * L)
x1 = np.arange(1, n, 1) # vector of values from 1 to 99
x2 = np.full(shape=n-1, fill_value = 1) # vector of 99 "1"
x3 = np.divide(x2,x1) # divied vector x2 by x1 (1/1, 1/2, 1/3 ...)
somme = np.sum(x3) # sum of the vector x3  
rho = np.divide(bin_breakpoints, somme) # get 4Ner * bin size 
r = np.divide(rho, 4 * (2.5*10**5) * 500) # get r in M/bp


### Store underlying and simulated recombination rates binned in 500 bp windows in a file
bin_edges=bin_edges[:-1].tolist()
y=y.tolist()
r=r.tolist()

rec = open('./input_LDhelmet/Underlying_Simulated_Map.txt', 'w')
for i in range(len(bin_edges)):
	pos = str(bin_edges[i]) # 1st column = positions
	real_rate = str(y[i]) # 2nd column = underlying rates
	simul_rate = str(r[i]) # 3rd column = simulated rates
	rec.write(pos + "\t" + real_rate + "\t" + simul_rate + "\n")

	
rec.close()

### Export tree sequence to a VCF file
with open("./input_LDhelmet/Admixture_pop_tmp.vcf", "w") as vcf_file:
    ts.write_vcf(vcf_file, ploidy=2) # precise the ploidie (= 2 for diploids) 

