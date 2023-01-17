#####################################################################
###############         HOTSPOT DETECTION         ###################
#####################################################################


########## DESCRIPTION ############
#### This script contains a function which called hotspots in a recombination map smoothed in 2.5 kb windows. 
## In the function, hotspots are defined as 2.5 kb windows (the resolution of the map) that are X time hotter than the 50kk flanking regions. 
## The threshold X to define the hotspots can be modified by the user (e.g. 2.5, 5, 10). 


########## HOTSTPOT FUNCTION ############
#### Function "hotspot" to call recombination hotspots
## It takes as argument a dataframe containing position, real, simul and infer rate, and the hotspot threshold value
hotspot <- function(rec_df, threshold) {
  AR <- rec_df
  HOT <- data.frame(matrix(NA, ncol=14))
  colnames(HOT) <- c("focal_win","wlo", "whi", "wmean_real", "wmean_simul", "wmean_infer", "rlo", "rhi", "rmean_real", "rmean_simul", "rmean_infer", "hot_real", "hot_simul", "hot_infer")
  for (i in 1:nrow(AR)){   # loop over the lines of the AR dataframe 
    HOT[i,1]<- i
    # set lower bound to a 2.5kb focal window centered on position "i" (used to calculate the recombination rate of the flanking region)
    wlo=(i-1)
    if (wlo<1)
    {wlo=1}
    HOT[i,2]<- wlo
    # set upper bound to a 2.5kb focal window centered on position "i" (used to calculate the recombination rate of the flanking region)
    whi=(i+1)
    if (whi>nrow(AR))
    {whi=nrow(AR)}
    HOT[i,3]<- whi
    # real recombination rate of position "i" 
    wmean_real <- AR[i,2]
    HOT[i,4]<- wmean_real
    # simul recombination rate of position "i" 
    wmean_simul <- AR[i,3]
    HOT[i,5]<- wmean_simul
    # infer recombination rate of position "i" 
    wmean_infer <- AR[i,4]
    HOT[i,6]<- wmean_infer
    # set lower bound to a flanking region of total length 100kb centered on position "i"
    rlo=(i-100000/(2*2500))
    if (rlo<1)
    {rlo=1}
    HOT[i,7]<- rlo
    # set upper bound to a flanking region of total length 100kb centered on position "i"
    rhi=(i+100000/(2*2500))
    if (rhi>nrow(AR))
    {rhi=nrow(AR)}
    HOT[i,8]<- rhi
    # mean recombination rate calculation for region of length 100kb centered on position "i", excluding positions in focal window for real_rate
    rmean_real=mean(na.omit(AR[c(rlo:wlo,whi:rhi),2]))
    HOT[i,9]<- rmean_real
    # mean recombination rate calculation for region of length 100kb centered on position "i", excluding positions in focal window for simul_rate
    rmean_simul=mean(na.omit(AR[c(rlo:wlo,whi:rhi),3]))
    HOT[i,10]<- rmean_simul
    # mean recombination rate calculation for region of length 100kb centered on position "i", excluding positions in focal window for infer_rate
    rmean_infer=mean(na.omit(AR[c(rlo:wlo,whi:rhi),4]))
    HOT[i,11]<- rmean_infer
    # output a vector of 3 TRUE/FALSE elements specifying hotspot status for Underlying/Simulated/Inferred recombination landscapes
    HOT[i,12] <- (wmean_real>threshold*rmean_real) 
    HOT[i,13] <- (wmean_simul>threshold*rmean_simul)
    HOT[i,14] <- (wmean_infer>threshold*rmean_infer)
    
  }
  HOT <- cbind(rec_df$positions,HOT)
  colnames(HOT)[1] <-"positions"
  return(HOT)
}




