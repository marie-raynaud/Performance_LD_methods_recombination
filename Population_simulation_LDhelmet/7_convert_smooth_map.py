########################################################################
################ CONVERT and SMOOTH RECOMBINATION RATES ################
########################################################################

########## DESCRIPTION ##########
### This script convert the recombination rates of the underlying and the simulated landscapes in units of rho by bp (scaling by Ne), and smooth the recombination maps in bins of 500 bp and 2.5 kb.
## It outputs two files: one file containing the bin positions, the underlying, simulated and inferred rates binned in 500 bp windows, and a second one binned in 2.5kb windows. 

########## VARIABLES ##########
## The bin size are hard-coded (500bp and 2.5kb) but can be modified by changing the number of bins


########## LOAD PACKAGES ##########
import scipy.stats
import numpy as np

### Load raw recombination map inferred by LDhelmet
data = np.loadtxt('~/Documents/Simulation_LDhelmet/output_LDhelmet/LDhelmet_raw_rates.txt')
positions = data[:,0] # retrieve positions of the left SNPs in a vector
rates = data[:,1] # retrieve recombination rates in a vector

### Smooth recombination rates in windows of 500 bp ###

num_bins = 2000 # define number of bins (2000 bins to generate 500 bp windows from a chromosomal segment of 1 Mb)
# Average recombination rates (values between adjacent SNPs) in 500 bp windows
bin_rates_500, bin_edges_500, _ = scipy.stats.binned_statistic(
    positions, rates, bins=num_bins, range=(0,1000000))

### Create an output file containing the underlying, simulated and inferred rates smoothed in 500 bp windows
# Load the file containing the underlying and simualted rates smoothed in 500 bp windows 
u_s_map = np.loadtxt('~/Documents/Simulation_LDhelmet/input_LDhelmet/Underlying_Simulated_Map.txt')
x1 = recap[:,0] # retrieve the positions
x2 = recap[:,1] # retrieve the underlying rates
x3 = recap[:,2] # retrieve the simualted rates 

### Convert the underlying and simulated rates, that are in units of M/bp, into units of rho/bp, by scaling by the effective size: rho=4Ne*r
x2 = x2 * (4 * (2.5*10**4)) # convert underlying rates into rho/bp 
x3 = x3 * (4 * (2.5*10**4)) # convert simulated rates into rho/bp

# Convert np.array into list
x1=x1.tolist() # position
x2=x2.tolist() # underlying rates
x3=x3.tolist()  # simulated rates
bin_rates_500=bin_rates_500.tolist() # inferred rates
bin_edges_500=bin_edges_500[:-1].tolist() # bin positions


### Store the bin positions, underlying, simulated and inferred rates in units of rho/bp smoothed in 500 bp windows in a text file
rec = open('Constant_size_pop_MAP_W500.txt', 'w')
for i in range(len(bin_edges_500)):
	pos = str(bin_edges_500[i]) # 1st column: bin positions
	real_rate = str(x2[i]) # 2nd column : underlying rates
	simul_rate = str(x3[i]) # 3rd column : simulated rates
	infer_rate = str(bin_rates_500[i]) # 4th column : inferred rates
	rec.write(pos + "\t" + real_rate + "\t" + simul_rate + "\t" + infer_rate +"\n")


rec.close()

### Smooth recombination rates in windows of 2500 bp ### 

num_bins = 400 # define number of bins (400 bins to generate 2500 bp windows from a chromosomal segment of 1 Mb)

# Average inferred recombination rates (values between adjacent SNPs) in 2500 bp windows
bin_rates_2500, bin_edges_2500, _ = scipy.stats.binned_statistic(
    positions, rates, bins=num_bins, range=(0,1000000))
	
# Average underlying recombination rates (converted in rho/bp) in 2500 bp windows
x2_2500, bin_edges_2500, _ = scipy.stats.binned_statistic(
    x1, x2, bins=num_bins, range=(0,1000000))

# Average simulated recombination rates (converted in rho/bp) in 2500 bp windows
x3_2500, bin_edges_2500, _ = scipy.stats.binned_statistic(
	x1, x3, bins=num_bins, range=(0,1000000))

# Convert np.array into list
x2_2500=x2_2500.tolist() # underlying rates
x3_2500=x3_2500.tolist() # simulated rates
bin_rates_2500=bin_rates_2500.tolist() # inferred rates
bin_edges_2500=bin_edges_2500[:-1].tolist() # bin positions


### Store the bin positions, underlying, simulated and inferred rates in units of rho/bp smoothed in 2500 bp windows in a text file
rec = open('Constant_size_pop_MAP_W2500.txt', 'w')
for i in range(len(bin_edges_2500)):
	pos = str(bin_edges_2500[i]) # 1st column: bin positions
	real_rate = str(x2_2500[i]) # 2nd column: underlying rates
	simul_rate = str(x3_2500[i]) # 3rd column: simulated rates
	infer_rate = str(bin_rates_2500[i]) # 4th column: inferred rates
	rec.write(pos + "\t" + real_rate + "\t" + simul_rate + "\t" + infer_rate +"\n")


rec.close()

quit()




