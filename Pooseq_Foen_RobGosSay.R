####  Analysis of PoolSeq data for Weber and Steinel ms, speed up with data.table
#### FP


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
library(data.table)

#################################################################
#  2.  Load data and organize as needed
#################################################################

# Load data

setwd("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/Data")

# These Fst values are corrected versions using "New_fst_GosRobSay_pools.R"

### Load Fst data
{
  RobSayFst <- fread("./allele_freq_Fst/RobSayFst_poolfstat.csv", header = F, sep = "\t", nrows = 6637336)
  GosRobFst <- fread("./allele_freq_Fst/GosRobFst_poolfstat.csv", header = F, sep = "\t", nrows = 6637336)
  GosSayFst <- fread("./allele_freq_Fst/GosSayFst_poolfstat.csv", header = F, sep = "\t", nrows = 6637336)
  GosSNPstats <- fread("./allele_freq_Fst/NewGosSNP_Stats.csv", header = T, nrows = 6637336)
  RobSNPstats <- fread("./allele_freq_Fst/NewRobSNP_Stats.csv", header = T, nrows = 6637336)
  SaySNPstats <- fread("./allele_freq_Fst/NewSaySNP_Stats.csv", header = T, nrows = 6637336)
  SaySNPstats <- SaySNPstats[,-1]
}


{
  #head(RobSayFst)
  #head(RobSNPstats)
  #head(SaySNPstats)
  
  FstAll <- data.table(GosSNPstats[,c(2:7)], RobSNPstats[,6:7], SaySNPstats[,6:7], GosRobFst[,3], GosSayFst[,3], RobSayFst[,3])
  head(FstAll)
  colnames(FstAll) <- c("LG", "Pos", "Base_A", "Base_a", "GosN_A", "GosNreads", "RobN_A", "RobNreads", "SayN_A", "SayNreads", "GRFst", "GSFst" , "RSFst")
  FstAll[,GosFreqA := FstAll$GosN_A/FstAll$GosNreads]
  FstAll[,RobFreqA := FstAll$RobN_A/FstAll$RobNreads]
  FstAll[,SayFreqA := FstAll$SayN_A/FstAll$SayNreads]
 
  # Clean up some memory
  rm(GosRobFst, GosSayFst, GosSNPstats, RobSayFst, RobSNPstats, SaySNPstats)
  
  ### Fst shouldn't be zero for the Rar calculations below            
  FstAll[,GRFst.nonneg := replace(FstAll$GRFst, which(FstAll$GRFst <= 0), NaN)]
  FstAll[,GSFst.nonneg := replace(FstAll$GSFst, which(FstAll$GSFst <= 0), NaN)]
  FstAll[,RSFst.nonneg := replace(FstAll$RSFst, which(FstAll$RSFst <= 0), NaN)]
  
  ### The Fst = Nan all reflect instances where allele frequencies are exactly the same in two pops so Fst should be = 0, to be conservative I set them to 0.001 rather than zero
  FstAll[is.nan(GRFst.nonneg), GRFst.nonneg := 0.00001]
  FstAll[is.nan(GSFst.nonneg), GSFst.nonneg := 0.00001]
  FstAll[is.nan(RSFst.nonneg), RSFst.nonneg := 0.00001]

  # When Fst = 1, when we do log(1-Fst), it becomes infinite so I replace 1.0 with 0.999  
  FstAll[GRFst.nonneg == 1.0, GRFst.nonneg := 0.99999]
  FstAll[GSFst.nonneg == 1.0, GSFst.nonneg := 0.99999]
  FstAll[RSFst.nonneg == 1.0, RSFst.nonneg := 0.99999]

  # Set a copy number threshold as 20 for the minimum and 150 for the maximum.
  min_cov<- 40
  max_cov<-150
  FstAll[,NreadsOK:=F]
  FstAll[GosNreads >=min_cov & RobNreads >=min_cov & SayNreads >=min_cov & GosNreads <=max_cov & RobNreads <=max_cov & SayNreads <=max_cov,
         NreadsOK:= T]
} 

#################################################################
  #  3.  Combined genome-wide plots of PoolSeq, Pi, Tajima's D
#################################################################
  
  # Caculate PBS for Rob branch, Gos branch, and Marine branch
  # From Sequencing of Fifty Human Exomes Reveals Adaptation to High Altitude Supplement
  FstAll[,pbs.r := (-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$RSFst.nonneg )) - (- log(1-FstAll$GSFst.nonneg ))/2]
  FstAll[,pbs.g := (-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$RSFst.nonneg ))/2]
  FstAll[,pbs.m := (-log(1-FstAll$RSFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$GRFst.nonneg ))/2]

  ## ********** NEED TO GET NUMERIC LG and cumulative position data for plotting 
  FstAll[substr(LG, 1, 2) == "MT", LGn := 0]
  FstAll[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  FstAll[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))] 
  
  FstAll[,PBS.percentile.unsmoothed:=ecdf(pbs.r)(pbs.r)]
  FstAll_top5<-FstAll[PBS.percentile.unsmoothed > 0.95,]
# Smoothed Fst and PBS
rm(list=setdiff(ls(), "FstAll"))

windowsize <- 1000 

# Put every SNP into a bin, bin 1000 include all the SNPs from 1-1000
FstAll[,SNP_window:= (Pos%/%windowsize+1)*windowsize]

# Smooth SNPs, only use the SNPs, which as reads number satisfy conditions
temp_BinSNPs<-FstAll[NreadsOK==T,.("PBS.r" = mean(pbs.r), "nSNPs" = .N),by=.(LGn, SNP_window)]

# Create a data.table with the full length of genome for temp_BinSNPs to join in.
chrlengths <- read.table("Stickle_chr_lengths.txt", sep = ",")
chr <- FstAll[, .(chr_length = max(Pos)), by = LGn]
chr[2:22,chr_length:= as.vector(as.matrix(chrlengths) )]
chr[,cumulative_chrlengths := cumsum(chr_length)]
chr[,cumulative_chrSTART := shift(cumulative_chrlengths,1,0,"lag")]
chr_forjoin<-chr[rep(1:.N,chr_length%/%windowsize+1)][,window_pos := ((1:.N))*windowsize, by = LGn]

# Join the smoothed Fst table with the chr_forjoin, the latter has the full length of chromosome
setkey(chr_forjoin, LGn, window_pos)
setkey(temp_BinSNPs, LGn, SNP_window)

PBS_binned <- temp_BinSNPs[chr_forjoin]
PBS_binned[,PBS.percentile.score:=1-ecdf(PBS.r)(PBS.r)]
PBS_binned[,cumulative_window_pos:=SNP_window+cumulative_chrSTART]
PBS_binned <- PBS_binned[,!c("chr_length","cumulative_chrlengths","cumulative_chrSTART")]
#BinSNPs[is.na(nSNPs),nSNPs:=0]
#BinSNPs[is.na(PBS.r),PBS.r:=0]

######################
#  Load Tajima's D 
TajD_r_filelist = list.files(path = "./TajD/subsample50x_mincount2_finished", pattern="Rob.Q15", full.names = T) 
TajD_r <- do.call(rbind,lapply(TajD_r_filelist,function(i){fread(i, sep = "\t", header = F)}))
colnames(TajD_r) <- c("LG", "Pos", "D.nSNPs", "x.r", "D.r")

TajD_r[substr(LG, 1, 2) == "MT", LGn := 0]
TajD_r[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
TajD_r[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
TajD_r<-TajD_r[,D.r := as.numeric(D.r)]
TajD_r[,D.percentile.score:=ecdf(D.r)(D.r)]
#TajD_r<-TajD_r[!is.na(D.r)]
setorder(TajD_r, LGn,Pos)

setkey(TajD_r, LGn, Pos)
setkey(PBS_binned, LGn, SNP_window)

TajD_PBS <- TajD_r[PBS_binned]
TajD_PBS_cleaned <- TajD_PBS[,.(LGn,Pos,cumulative_window_pos,D.nSNPs,D.r, D.percentile.score,"PBS.nSNPs" = nSNPs, PBS.r, PBS.percentile.score)]

#############
# Load Pi
Pi_r_filelist = list.files(path = "./per_chromosom_pi", pattern="Rob.Q15", full.names = T) 
Pi_r <- do.call(rbind,lapply(Pi_r_filelist,function(i){fread(i, sep = "\t", header = F)}))
colnames(Pi_r) <- c("LG", "Pos", "Pi.nSNPs", "x.r", "Pi.r")

Pi_r[substr(LG, 1, 2) == "MT", LGn := 0]
Pi_r[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
Pi_r[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
Pi_r<-Pi_r[,Pi.r := as.numeric(Pi.r)]
Pi_r[,Pi.percentile.score:=ecdf(Pi.r)(Pi.r)]
#Pi_r<-Pi_r[!is.na(Pi.r)]
setorder(Pi_r, LGn,Pos)

setkey(Pi_r, LGn, Pos)
setkey(TajD_PBS_cleaned, LGn, Pos)

# rolljoin in data.table can assign values to its closest neighbor. -Inf means backward rolling. Note that Pi.nSNPs is the sum of those 2-3 rows
Pi_TajD_PBS <- Pi_r[TajD_PBS_cleaned,roll = -Inf]
Pi_TajD_PBS[,combined.percentile.score:=Pi.percentile.score*D.percentile.score*PBS.percentile.score]
Pi_TajD_PBS[,combined.percentile:=ecdf(combined.percentile.score)(combined.percentile.score)]
setorder(Pi_TajD_PBS,combined.percentile)
Pi_TajD_PBS_cleaned <- Pi_TajD_PBS[,.(LGn,Pos,cumulative_window_pos, 
                                      Pi.nSNPs, Pi.r, Pi.percentile.score,
                                      D.nSNPs, D.r, D.percentile.score,
                                      PBS.nSNPs, PBS.r, PBS.percentile.score,
                                      combined.percentile.score,combined.percentile)]
combined_top5<-Pi_TajD_PBS_cleaned[combined.percentile<0.05,]

fwrite(combined_top5, "Pi_TajD_PBS_smoothed_top5percent.csv")
fwrite(Pi_TajD_PBS_cleaned, "Pi_TajD_PBS_smoothed_full.csv")


#################################################################
#  4.  Identify main sites of interest
#################################################################

source("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/R Code/stickleback_poolseq_RobGosSay/gwscaR_plot.R")
combined_top5_df <- as.data.frame(combined_top5)
FstAll_top5_df <- as.data.frame(FstAll_top5[FstAll_top5$NreadsOK==T,])
Pi_TajD_PBS_cleaned_df <- as.data.frame(Pi_TajD_PBS_cleaned)

fst.plot(fst.dat = FstAll_top5_df ,scaffold.widths=chr[,.(LGn,chr_length)],scaffs.to.plot= 1:21,
        fst.name="pbs.r", chrom.name="LGn", bp.name="Pos",
        y.lim=c(8,25),xlabels=1:21,xlab.indices=NULL,
        axis.size=0.5,pt.cols=c("darkgrey","lightgrey"),pt.cex=0.5)
