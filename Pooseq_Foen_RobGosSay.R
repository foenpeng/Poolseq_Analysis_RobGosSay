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
#library(ape)

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


#####
# Merge Fst info and calculate PBS and Kirkpatrick's Rar:
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

  # Set a copy number threshold as 20 the minimum and 150 the maximum.
  FstAll[,NreadsOK:=F]
  FstAll[GosNreads >19 & RobNreads >19 & SayNreads >19 & GosNreads <151 & RobNreads <151 & SayNreads <151,
         NreadsOK:= T]
  
  # Caculate PBS for Rob branch, Gos branch, and Marine branch
  # From Sequencing of Fifty Human Exomes Reveals Adaptation to High Altitude Supplement
  FstAll[,pbs.r := (-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$RSFst.nonneg )) - (- log(1-FstAll$GSFst.nonneg ))/2]
  FstAll[,pbs.g := (-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$RSFst.nonneg ))/2]
  FstAll[,pbs.m := (-log(1-FstAll$RSFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$GRFst.nonneg ))/2]

  #########
  ### ********** NEED TO GET NUMERIC LG and cumulative position data for plotting 
  FstAll[substr(LG, 1, 2) == "MT", LGn := 0]
  FstAll[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  FstAll[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))] 
  
  #### Define a continuous genomic x axis, MTDNA near zero
  Chrlengths <- read.table("Stickle_chr_lengths.txt", sep = ",")
  chrmax <- FstAll[, .(max_length = max(Pos)), by = LGn]
  chrmax[2:22,max_length:= as.vector(as.matrix(Chrlengths) )]
  chrmax[,Cumulative_Chrlengths := cumsum(max_length)]
  Cumulative_ChrSTART <- chrmax[,shift(Cumulative_Chrlengths,1,0,"lag")]
  FstAll[, Cumulative_position := Pos + Cumulative_ChrSTART[FstAll[,LGn]]]
  

}

##############
# Smoothed Fst and PBS
rm(list=setdiff(ls(), "FstAll"))
windowsize <- 1000 

FstAll[,SNP_window:= (Pos%/%windowsize)*windowsize]
temp_BinSNPs<-FstAll[NreadsOK==T,.("PBS.r" = mean(pbs.r), "nSNPs" = .N),by=.(LGn, SNP_window)]

# Create a data.table with the full length of genome for BinSNPs to join in.
Chrlengths <- read.table("Stickle_chr_lengths.txt", sep = ",")
chr <- FstAll[, .(chr_length = max(Pos)), by = LGn]
chr[2:22,chr_length:= as.vector(as.matrix(Chrlengths) )]
chr[,Cumulative_Chrlengths := cumsum(chr_length)]
chr[,Cumulative_ChrSTART := shift(Cumulative_Chrlengths,1,0,"lag")]
chr_forjoin<-chr[rep(1:.N,chr_length%/%windowsize+1)][,window_pos := ((1:.N)-1)*windowsize, by = LGn]

setkey(chr_forjoin, LGn, window_pos)
setkey(temp_BinSNPs, LGn, SNP_window)

BinSNPs <- temp_BinSNPs[chr_forjoin]
BinSNPs[,Cumulative_window_pos:=SNP_window+Cumulative_ChrSTART]
BinSNPs <- BinSNPs[,!c("chr_length","Cumulative_Chrlengths","Cumulative_ChrSTART")]
BinSNPs[is.na(nSNPs),nSNPs:=0]
BinSNPs[is.na(PBS.r),PBS.r:=0]

######################
#  Load Tajima's D 
setwd("./TajD/subsample50x_mincount2_finished/")
TajD_filelist = list.files( pattern="group.*Rob.Q15") 
TajD.r <- do.call(rbind,lapply(TajD_filelist,function(i){fread(i, , sep = "\t", header = F)}))
colnames(TajD.r) <- c("LG", "Pos", "Nsnp.r", "x.r", "D.r")

TajD.r[substr(LG, 1, 2) == "MT", LGn := 0]
TajD.r[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
TajD,r[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))] 
