####  July 30 2018 Analysis of PoolSeq data for Weber and Steinel ms
#### DIB


### PURPOSE:
# 1. Establish that selection has acted on sites in Roberts, Gosling, and Marine populations.
# 2. Identify where those sites are, and what genes are in them

#################################################################
##### CONTENTS:
# 1. Load packages
# 2. Load data
# 3. Combined genome-wide plots of PoolSeq, Pi, Tajima's D
# 4. Identify main sites of interest
# 5. Zoomed-in on sites of interest
#################################################################


#################################################################
#  NOTES:

#################################################################






#################################################################
#  1.  Load packages
#################################################################

library(ape)









#################################################################
#  2.  Load data and organize as needed
#################################################################

# Load data

setwd("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/Data")
# These Fst values are corrected versions using "New_fst_GosRobSay_pools.R"

### Load Fst data
{
RobSayFst <- read.csv("./allele_freq_Fst/RobSayFst_poolfstat.csv", header = F, sep = "\t", nrows = 6637336, colClasses = c("factor", "integer", "numeric"))
GosRobFst <- read.csv("./allele_freq_Fst/GosRobFst_poolfstat.csv", header = F, sep = "\t", nrows = 6637336, colClasses = c("factor", "integer", "numeric"))
GosSayFst <- read.csv("./allele_freq_Fst/GosSayFst_poolfstat.csv", header = F, sep = "\t", nrows = 6637336, colClasses = c("factor", "integer", "numeric"))
GosSNPstats <- read.csv("./allele_freq_Fst/NewGosSNP_Stats.csv", header = T, nrows = 6637336)
RobSNPstats <- read.csv("./allele_freq_Fst/NewRobSNP_Stats.csv", header = T, nrows = 6637336)
SaySNPstats <- read.csv("./allele_freq_Fst/NewSaySNP_Stats.csv", header = T, nrows = 6637336)
SaySNPstats <- SaySNPstats[,-1]
}


#####
# Merge Fst info and calculate PBS and Kirkpatrick's Rar:
{
head(RobSayFst)
head(RobSNPstats)
head(SaySNPstats)

FstAll <- data.frame(GosSNPstats[,c(2:7)], RobSNPstats[,6:7], SaySNPstats[,6:7], GosRobFst[,3], GosSayFst[,3], RobSayFst[,3])
head(FstAll)
colnames(FstAll) <- c("LG", "Pos", "Base_A", "Base_a", "GosN_A", "GosNreads", "RobN_A", "RobNreads", "SayN_A", "SayNreads", "GRFst", "GSFst" , "RSFst")
FstAll$GosFreqA <- FstAll$GosN_A/FstAll$GosNreads
FstAll$RobFreqA <- FstAll$RobN_A/FstAll$RobNreads
FstAll$SayFreqA <- FstAll$SayN_A/FstAll$SayNreads

### Fst shouldn't be zero for the Rar calculations below            ??? should convert negative Fst to zeros

FstAll$GRFst.nonneg <- replace(FstAll$GRFst, which(FstAll$GRFst < 0), NaN)
FstAll$GSFst.nonneg <- replace(FstAll$GSFst, which(FstAll$GSFst < 0), NaN)
FstAll$RSFst.nonneg <- replace(FstAll$RSFst, which(FstAll$RSFst < 0), NaN)

### The Fst = Nan all reflect instances where allele frequencies are exactly the same in two pops so Fst should be = 0, to be conservative I set them to 0.001 rather than zero
FstAll$GRFst.nonneg[is.nan(FstAll$GRFst.nonneg )] <- 0.001
FstAll$GSFst.nonneg[is.nan(FstAll$GSFst.nonneg )] <- 0.001
FstAll$RSFst.nonneg[is.nan(FstAll$RSFst.nonneg )] <- 0.001

# When Fst = 1, Kirkpatrick's Rar becomes NaN when we do log(1-Fst) so I replace 1.0 with 0.999
FstAll$GRFst.nonneg[FstAll$GRFst.nonneg == 1.0] <- 0.999
FstAll$GSFst.nonneg[FstAll$GSFst.nonneg == 1.0] <- 0.999
FstAll$RSFst.nonneg[FstAll$RSFst.nonneg == 1.0] <- 0.999

# Caculate PBS for Rob branch, Gos branch, and Marine branch
# From Sequencing of Fifty Human Exomes Reveals Adaptation to High Altitude Supplement
FstAll$pbs.r <- (-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$RSFst.nonneg )) - (- log(1-FstAll$GSFst.nonneg ))/2
FstAll$pbs.g <- (-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$RSFst.nonneg ))/2
FstAll$pbs.m <- (-log(1-FstAll$RSFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$GRFst.nonneg ))/2

### Now add kirkpatrick's Ra
FstAll$Ra.r <- (FstAll$GRFst.nonneg  + FstAll$RSFst.nonneg  - FstAll$GSFst.nonneg )/FstAll$GSFst.nonneg 
FstAll$Ra.g <-  (FstAll$GRFst.nonneg  + FstAll$GSFst.nonneg  - FstAll$RSFst.nonneg )/FstAll$RSFst.nonneg 
FstAll$Ra.m <-  (FstAll$GSFst.nonneg  + FstAll$RSFst.nonneg  - FstAll$GRFst.nonneg )/FstAll$GRFst.nonneg 

# Now remove NAs and Infs
FstAll$Ra.r[is.nan(FstAll$Ra.r)] <- NA
FstAll$Ra.g[is.nan(FstAll$Ra.g)] <- NA
FstAll$Ra.m[is.nan(FstAll$Ra.m)] <- NA

  #  inf vals involve Fst = # / Fst = 0  Not really NA, these are actually quite interesting. So note them, then change to NA for plotting without losing that info
FstAll$Ra.r_isINF <- is.infinite(FstAll$Ra.r)
FstAll$Ra.g_isINF <- is.infinite(FstAll$Ra.g)
FstAll$Ra.m_isINF <- is.infinite(FstAll$Ra.m)

FstAll$Ra.r[is.infinite(FstAll$Ra.r)] <- max(FstAll$Ra.r, na.m = T)
FstAll$Ra.g[is.infinite(FstAll$Ra.g)] <- max(FstAll$Ra.g, na.m = T)
FstAll$Ra.m[is.infinite(FstAll$Ra.m)] <- max(FstAll$Ra.m, na.m = T)

#########
### ********** NEED TO GET NUMERIC LG and cumulative position data for plotting     
FstAll$LG <- as.character(FstAll$LG)
FstAll$isMT <- substr(FstAll$LG, 1, 2) == "MT"
FstAll$isLG <- substr(FstAll$LG, 1, 2) == "gr"
FstAll$isScaffold <- substr(FstAll$LG, 1, 2) == "sc"
FstAll$LGn <- substr(FstAll$LG, 6, nchar(FstAll$LG))
FstAll$LGn[FstAll$isMT | FstAll$isScaffold] <- NA
FstAll$LGn <- as.numeric(as.roman(FstAll$LGn))
FstAll$LGn[FstAll$isMT] <- 0
FstAll$LGn[FstAll$isScaffold] <- as.numeric(substr(FstAll$LG, 10, nchar(FstAll$LG))[FstAll$isScaffold])


#### Define a continuous genomic x axis, MTDNA near zero
Chrlengths <- read.table("Stickle_chr_lengths.txt", sep = ",")
chrmax <- tapply(FstAll$Pos, FstAll$LGn, max)
plot(chrmax[2:22] ~ t(Chrlengths))
chrmax_names <- names(chrmax)
chrmax <- as.vector(chrmax)
chrmax[2:22] <- as.vector(as.matrix(Chrlengths) )
Cumulative_Chrlengths <- cumsum(chrmax)
names(Cumulative_Chrlengths) <- chrmax_names
Cumulative_ChrSTART <- c(0, Cumulative_Chrlengths[-length(Cumulative_Chrlengths)])
names(Cumulative_ChrSTART) <- chrmax_names

FstAll$Cumulative_position <- FstAll$Pos + Cumulative_ChrSTART[as.character(FstAll$LGn)]

}



{
  # Clean up some memory
  rm(GosRobFst, GosSayFst, GosSNPstats, RobSayFst, RobSNPstats, SaySNPstats)
}





##############
# Smoothed Fst and PBS and RaR
windowsize <- 10000 #
# create blank dataframe for saving results
FstsSmooth <- as.data.frame(matrix(ncol = 12, nrow = sum(chrmax)/(0.9*windowsize )  ))
  counter <- 0
  # Iterate through all chromosomes / scaffolds
for(chr in 1:length(chrmax)){
    # Check to be sure scaffold is long enough to us
  if(chrmax[chr] > windowsize){
      # define start of window
    startsitelist <- seq(0, chrmax[chr] - windowsize, by = windowsize)
    LGname <- chrmax_names[chr]
    print(LGname)
    for(startsite in startsitelist){
    counter <- counter + 1
      # acquire data from within focal window
    FstAll.focal <- FstAll[FstAll$LGn == as.numeric(LGname) & FstAll$Pos > startsite & FstAll$Pos < startsite + windowsize ,  ]
      # Add averaged data from window into dataframe
    FstsSmooth[counter,] <-data.frame(LG = as.numeric(LGname), 
                                      Pos = startsite + windowsize/2, 
                                      GRFst.nonneg = mean(FstAll.focal$GRFst.nonneg, na.rm = T),  
                                      GSFst.nonneg = mean(FstAll.focal$GSFst.nonneg, na.rm = T),  
                                      RSFst.nonneg = mean(FstAll.focal$RSFst.nonneg, na.rm = T), 
                                      PBS.r = mean(FstAll.focal$pbs.r, na.rm= T),
                                      PBS.g = mean(FstAll.focal$pbs.g, na.rm= T),
                                      PBS.m = mean(FstAll.focal$pbs.m, na.rm= T), 
                                      Ra.r = mean(FstAll.focal$Ra.r, na.rm = T),  
                                      Ra.g = mean(FstAll.focal$Ra.g, na.rm = T),    
                                      Ra.m = mean(FstAll.focal$Ra.m, na.rm = T), 
                                      nSNPs = nrow(FstAll.focal))
  } }   
}
  
colnames(FstsSmooth) <- c("LG", "Pos", "GRFst.nonneg", "GSFst.nonneg", "RSFst.nonneg", "PBS.r", "PBS.g", "PBS.m", "Ra.r", "Ra.g", "Ra.m", "nSNPs")
FstsSmooth
dim(FstsSmooth)
FstsSmooth.cleaned <- FstsSmooth[FstsSmooth$nSNPs >0,]
dim(FstsSmooth.cleaned)
#FstsSmooth.cleaned <- FstsSmooth.cleaned[rowSums(is.na(FstsSmooth.cleaned)) == 0, ]
# make culumative x values
FstsSmooth.cleaned$Cumulative_position <- FstsSmooth.cleaned$Pos + Cumulative_ChrSTART[as.character(FstsSmooth.cleaned$LG)]
# Save
write.csv(FstsSmooth.cleaned, "./allele_freq_Fst/Binned10000bpFst_and_Ra.csv")



# Exploratory plotting, can skip
# chrcol <- 2+2*as.numeric(FstsSmooth.cleaned$LG %%2 == 1)
# plot(Ra.r ~ Cumulative_position, FstsSmooth.cleaned, col = chrcol, type = "l", ylim = c(-10, 200000))
# plot(Ra.g ~ Cumulative_position, FstsSmooth.cleaned, col = chrcol, type = "l")
# plot(Ra.r ~ Cumulative_position, FstAll, col = chrcol, type = "l", ylim = c(-10, 10000000))
# 
# # Add in outlier points
# pointthreshold <- quantile(FstAll$Ra.r, 0.999)
# points(FstAll$Ra.r[FstAll$Ra.r > pointthreshold] ~ FstAll$Cumulative_position[FstAll$Ra.r > pointthreshold], FstAll, col = chrcol)






######################
#  Load Tajima's D 
{
{
TajD.g <- read.csv("./TajD/subsample50x_mincount2_finished/grouPi.ros.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F)
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupIII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupIV.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupV.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupVI.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupVII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupVIII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupIX.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupX.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXI.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXIII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXIV.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXV.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXVI.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXVII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXVIII.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXIX.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXX.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
TajD.g <- rbind(TajD.g, read.csv("./TajD/subsample50x_mincount2_finished/groupXXI.Gos.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
}

{
  TajD.r <- read.csv("./TajD/subsample50x_mincount2_finished/groupI.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F)
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupIII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupIV.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupV.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupVI.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupVII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupVIII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupIX.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupX.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXI.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXIII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXIV.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXV.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXVI.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXVII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXVIII.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXIX.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXX.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.r <- rbind(TajD.r, read.csv("./TajD/subsample50x_mincount2_finished/groupXXI.Rob.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
}

{
  TajD.m <- read.csv("./TajD/subsample50x_mincount2_finished/groupI.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD.txt", sep = "\t", header = F)
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupIII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupIV.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupV.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupVI.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupVII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupVIII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupIX.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupX.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXI.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXIII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXIV.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXV.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXVI.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXVII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXVIII.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXIX.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXX.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
  TajD.m <- rbind(TajD.m, read.csv("./TajD/subsample50x_mincount2_finished/groupXXI.Say.Q15.NoIndel.Broad.minCount2.minQ15.minCovFrac.5.2000win.1000step.tajD", sep = "\t", header = F))
}
  
  
  {
    TajD.g$site <- paste(TajD.g$V1, TajD.g$V2, sep = "_")
    TajD.r$site <- paste(TajD.r$V1, TajD.r$V2, sep = "_")
    TajD.m$site <- paste(TajD.m$V1, TajD.m$V2, sep = "_")
    
  }
  
{
  TajD <- merge(TajD.g, TajD.r, by ="site", all = T)
TajD <- merge(TajD, TajD.m, by ="site", all = T)
TajD <- TajD[,c(2:6, 9:11, 14:16)]
colnames(TajD) <- c("LG", "Pos", "Nsnp.g", "x.g", "D.g", "Nsnp.r", "x.r", "D.r", "Nsnp.m", "x.m", "D.m")
}

  {  
#########
### ********** NEED TO GET NUMERIC LG and cumulative position data for plotting
TajD$LG <- as.character(TajD$LG)
TajD$isMT <- substr(TajD$LG, 1, 2) == "MT"
TajD$isLG <- substr(TajD$LG, 1, 2) == "gr"
TajD$isScaffold <- substr(TajD$LG, 1, 2) == "sc"
TajD$LGn <- substr(TajD$LG, 6, nchar(TajD$LG))
TajD$LGn[TajD$isMT | TajD$isScaffold] <- NA
TajD$LGn <- as.numeric(as.roman(TajD$LGn))
TajD$LGn[TajD$isMT] <- 0
#TajD$LGn[TajD$isScaffold] <- as.numeric(substr(TajD$LG, 10, nchar(TajD$LG))[TajD$isScaffold])


#### Define a continuous genomic x axis, MTDNA near zero
TajD$Cumulative_position <- TajD$Pos + Cumulative_ChrSTART[as.character(TajD$LGn)]
}

write.csv(TajD, "./TajD/MergedTajD.csv")
}
  
  
#######################
#   Load Pi
{
{
Pi.g <- read.csv("./per_chromosom_pi/groupI.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi.csv", sep = "\t", header = F)
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupIII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupIV.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupV.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupVI.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupVII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupVIII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupIX.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupX.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXI.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXIII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXIV.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXV.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXVI.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXVII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXVIII.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXIX.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXX.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
Pi.g <- rbind(Pi.g, read.csv("./per_chromosom_pi/groupXXI.Gos.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
}



{
  Pi.m <- read.csv("./per_chromosom_pi/groupI.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F)
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupIII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupIV.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupV.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupVI.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupVII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupVIII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupIX.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupX.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXI.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXIII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXIV.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXV.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXVI.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXVII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXVIII.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXIX.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXX.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.m <- rbind(Pi.m, read.csv("./per_chromosom_pi/groupXXI.Say.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
}


{
  Pi.r <- read.csv("./per_chromosom_pi/groupI.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F)
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupIII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupIV.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupV.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupVI.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupVII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupVIII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupIX.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupX.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXI.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXIII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXIV.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXV.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXVI.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXVII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXVIII.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXIX.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXX.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
  Pi.r <- rbind(Pi.r, read.csv("./per_chromosom_pi/groupXXI.Rob.Q15.NoIndel.Broad.minCount1.minCov40.5000win.2500step.pi", sep = "\t", header = F))
}

  {
    Pi.g$site <- paste(Pi.g$V1, Pi.g$V2, sep = "_")
    Pi.r$site <- paste(Pi.r$V1, Pi.r$V2, sep = "_")
    Pi.m$site <- paste(Pi.m$V1, Pi.m$V2, sep = "_")
    
  }
  
  {
    Pi <- merge(Pi.g, Pi.r, by ="site", all = T)
    Pi <- merge(Pi, Pi.m, by ="site", all = T)
    Pi <- Pi[,c(2:6, 9:11, 14:16)]
    colnames(Pi) <- c("LG", "Pos", "Nsnp.g", "x.g", "Pi.g", "Nsnp.r", "x.r", "Pi.r", "Nsnp.m", "x.m", "Pi.m")
  }
  
  {  
    #########
    ### ********** NEED TO GET NUMERIC LG and cumulative position data for plotting
    Pi$LG <- as.character(Pi$LG)
    Pi$isMT <- substr(Pi$LG, 1, 2) == "MT"
    Pi$isLG <- substr(Pi$LG, 1, 2) == "gr"
    Pi$isScaffold <- substr(Pi$LG, 1, 2) == "sc"
    Pi$LGn <- substr(Pi$LG, 6, nchar(Pi$LG))
    Pi$LGn[Pi$isMT | Pi$isScaffold] <- NA
    Pi$LGn <- as.numeric(as.roman(Pi$LGn))
    Pi$LGn[Pi$isMT] <- 0
    #Pi$LGn[Pi$isScaffold] <- as.numeric(substr(Pi$LG, 10, nchar(Pi$LG))[Pi$isScaffold])
    
    
    #### Define a continuous genomic x axis, MTDNA near zero
    Pi$Cumulative_position <- Pi$Pos + Cumulative_ChrSTART[as.character(Pi$LGn)]
  }
  
  write.csv(Pi, "./per_chromosom_pi/MergedPi.csv")
}


#################################################################
#  3.  Combined genome-wide plots of PoolSeq, Pi, Tajima's D
#################################################################


# 3. Combined genome-wide plots of PoolSeq, Pi, Tajima's D
# 4. Identify main sites of interest
# 5. Zoomed-in on sites of interest

















## Plot Fst values, Ra values, depth of coverage # Smoothing, but how to handle varying dpth of coverag if smoothed














ylimRob <- 300
ylimGos <- 450
ylimSay <- 600
{
{
par(mar = c(2,5,0.5,1), mfrow = c(3,1), oma = c(3,0,0,0))
{
plot(Ra.r ~ I(start/1000000), fsts.m, xlab = "Gene position (mb)", ylab = "PBS (Roberts Lake)", cex.lab = 1.5, pch = 16, ylim = c(0,ylimRob))
focalfsts <- fsts.m[fsts.m$Ra.r > 35,]
text(x = I(focalfsts$start/1000000), y = focalfsts$Ra.r, focalfsts$Associated.Gene.Name, pos = 4,cex = 1.5 )
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3675233/
text(x = 1761527/1000000, y = 45.56958, "Asap1/2", pos = 4,cex = 1.5 )

temp <- na.omit(fsts.m)
spl <- smooth.spline(x = temp$start/1000000, y = temp$Ra.r, spar = 0.3 )
lines(spl, lwd = 3)
}

{
  plot(Ra.g ~ I(start/1000000), fsts.m, xlab = "Gene position (mb)", ylab = "PBS (Gosling Lake)", cex.lab = 1.5, pch = 16, ylim = c(0,ylimGos))
  focalfsts <- fsts.m[fsts.m$Ra.g > 35,]
  text(x = I(focalfsts$start/1000000), y = focalfsts$Ra.g, focalfsts$Associated.Gene.Name, pos = 4,cex = 1.5 )
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3675233/
  temp <- na.omit(fsts.m)
  spl <- smooth.spline(x = temp$start/1000000, y = temp$Ra.g, spar = 0.3 )
  lines(spl, lwd = 3)
}

{
  plot(Ra.m ~ I(start/1000000), fsts.m, xlab = "Gene position (mb)", ylab = "PBS (Sayward)", cex.lab = 1.5, pch = 16, ylim = c(0,ylimSay))
  focalfsts <- fsts.m[fsts.m$Ra.m > 25,]
  text(x = I(focalfsts$start/1000000), y = focalfsts$Ra.m, focalfsts$Associated.Gene.Name, pos = 4,cex = 1.5 )
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3675233/
  temp <- na.omit(fsts.m)
  spl <- smooth.spline(x = temp$start/1000000, y = temp$Ra.m, spar = 0.3 )
  lines(spl, lwd = 3)
}
  mtext( "Gene chromosomal position, Mb", side = 1, outer = T, line = 1, cex = 1.5)
}

  
  
fsts.m[order(fsts.m$Ra.r, decreasing = T),][1:10,]
fsts.m <- fsts.m[-1,]
{
distmatrix <- fsts.m[fsts.m$Associated.Gene.Name == "mafbb",6:8][1,]
distmatrix <- c(0,distmatrix[1],distmatrix[3],distmatrix[1],0,distmatrix[2],distmatrix[3],distmatrix[2],0)
distmatrix <- matrix(as.numeric(distmatrix), nrow = 3, byrow = T)
rownames(distmatrix) = colnames(distmatrix) <- c("Roberts", "Gosling", "Marine")
distmatrix <- as.dist(distmatrix)
disttree <- nj(distmatrix)
distree.r <- root(disttree, 3)
par(fig = c(0.78, 0.98, 0.85, 0.99), new = T)
plot(distree.r, type = "u", cex = 0.01, label.offset = 3, edge.width = 2, no.margin = T, frame = T, rotate.tree = 30)
tiplabels(pch = c(16,17,1), cex = 2.2, adj = c(0.5, 0.5))
add.scale.bar(length = 0.1)
mtext("mafbb", side = 1, line = -1)

}


{
  distmatrix <- fsts.m[fsts.m$Associated.Gene.Name == "golga7",6:8][1,]
  distmatrix <- c(0,distmatrix[1],distmatrix[3],distmatrix[1],0,distmatrix[2],distmatrix[3],distmatrix[2],0)
  distmatrix <- matrix(as.numeric(distmatrix), nrow = 3, byrow = T)
  distmatrix
  rownames(distmatrix) = colnames(distmatrix) <- c("Roberts", "Gosling", "Marine")
  distmatrix <- as.dist(distmatrix)
      #distmatrix <- log(distmatrix)
  disttree <- nj(distmatrix)
  distree.r <- root(disttree, 3)
  par(fig = c(0.78,  0.98, 0.18, 0.32), new = T)
  plot(distree.r, type = "u", cex = 0.01, label.offset = 3, edge.width = 2, no.margin = T, frame = T, rotate.tree = 90)
  tiplabels(pch = c(16,17,1), cex = 2.2, adj = c(0.5, 0.50))
  mtext("golga7", side = 1, line = -1.2)
}

}


# Combine Fst, Pi and D ---------------------------------------------------

# combine value 
stepmean <- function(value,step){
  n<-length(value)
  remainder = n%%step
  result<-colMeans(matrix(value[1:(n-remainder)],step))
  result<-c(result,mean(value[(n-remainder+1):n]))
  return(result)
}