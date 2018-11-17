####  Analysis of PoolSeq data for Weber and Steinel ms, speed up with data.table
#### FP

# To-do
# 10.27, Chi-square test, check SNP difference among statistics, graph, annotation
# it seems necessary to set threshold value for pbs.m and pbs.g, becasue if those two are really high, pbs.r tend to be high
# chi-square test filtered those in top 5%

### PURPOSE:
# 1. Establish that selection has acted on sites in Roberts, Gosling, and Marine populations.
# 2. Identify where those sites are, and what genes are in them

#################################################################
##### CONTENTS:
# 1. Load packages
# 2. Load data
# 3. Calculate PBS and smoothing
# 4. Combined genome-wide plots of PoolSeq, Pi, Tajima's D
# 5. Identify main sites of interest
# 6. Zoomed-in on sites of interest
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

### Pre-proceesing Fst data
{
  FstAll <- data.table(GosSNPstats[,c(2:7)], RobSNPstats[,6:7], SaySNPstats[,6:7], GosRobFst[,3], GosSayFst[,3], RobSayFst[,3])
  head(FstAll)
  colnames(FstAll) <- c("LG", "Pos", "Base_A", "Base_a", "GosN_A", "GosNreads", "RobN_A", "RobNreads", "SayN_A", "SayNreads", "GRFst", "GSFst" , "RSFst")
  FstAll[,GosFreqA := FstAll$GosN_A/FstAll$GosNreads]
  FstAll[,RobFreqA := FstAll$RobN_A/FstAll$RobNreads]
  FstAll[,SayFreqA := FstAll$SayN_A/FstAll$SayNreads]

  # Clean up some memory
  rm(GosRobFst, GosSayFst, GosSNPstats, RobSayFst, RobSNPstats, SaySNPstats)
  
  # The Fst = Nan all reflect instances where allele frequencies are exactly the same in two pops so Fst should be = 0, 
  # Set all the NA and negative Fst as 0.0001 ~6000 (0.1% of all) SNPs of GSFst have value between 0 and 0.0001.
  FstAll[,GRFst.nonneg := replace(GRFst, which(GRFst < 0 | is.na(GRFst)), 0)]
  FstAll[,GSFst.nonneg := replace(GSFst, which(GSFst < 0 | is.na(GSFst)), 0)]
  FstAll[,RSFst.nonneg := replace(RSFst, which(RSFst < 0 | is.na(RSFst)), 0)]

  # When Fst = 1, when we do log(1-Fst), it becomes infinite so I replace 1.0 with 0.99
  # If I set it to 0.999 or 0.9999, it will create big PBS, because -log(1-0.9999) is a big value
  FstAll[GRFst.nonneg == 1.0, GRFst.nonneg := 0.99]
  FstAll[GSFst.nonneg == 1.0, GSFst.nonneg := 0.99]
  FstAll[RSFst.nonneg == 1.0, RSFst.nonneg := 0.99]

  # Set a copy number threshold as 20 for the minimum and 200 for the maximum. And mark them.
  # if I set max at 150, there will be a lot SNPs which pass Chi square test but filtered because of read numebrs.
  min_cov<- 20
  max_cov<-200
  FstAll[,NreadsOK:=F]
  FstAll[GosNreads >=min_cov & RobNreads >=min_cov & SayNreads >=min_cov & GosNreads <=max_cov & RobNreads <=max_cov & SayNreads <=max_cov,
         NreadsOK:= T]
} 

#################################################################
  #  3.  Calculate PBS and smoothing
#################################################################
  
  # Caculate PBS for Rob branch, Gos branch, and Marine branch
  # From Sequencing of Fifty Human Exomes Reveals Adaptation to High Altitude Supplement
  FstAll[,pbs.r := ((-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$RSFst.nonneg )) - (- log(1-FstAll$GSFst.nonneg )))/2]
  FstAll[,pbs.g := ((-log(1-FstAll$GRFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$RSFst.nonneg )))/2]
  #FstAll[,pbs.s := ((-log(1-FstAll$RSFst.nonneg )) + (- log(1-FstAll$GSFst.nonneg )) - (- log(1-FstAll$GRFst.nonneg )))/2]
  FstAll[,pbs.r.percentile:=ecdf(pbs.r)(pbs.r)]
  
  # Ra calculation
  FstAll[,Ra.r := (GRFst.nonneg  + RSFst.nonneg  - GSFst.nonneg)/GSFst.nonneg ]
  FstAll[,Ra.r.percentile:=ecdf(Ra.r)(Ra.r)]
  
  # Calculate dxy between Gos and Rob fishes.dxy = (p1∗(1−p2))+(p2∗(1−p1)). According to "Comparative analysis examining patterns of genomic differentiation across multiple episodes of population divergence in birds" 
  # FstAll[,dxy.GR := GosFreqA  * (1- RobFreqA)  + RobFreqA * (1- GosFreqA)] 
  # FstAll[,dxy.GR.percentile:=ecdf(dxy.GR)(dxy.GR)]
  
  # GET NUMERIC LG
  FstAll[substr(LG, 1, 2) == "MT", LGn := 0]
  FstAll[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  FstAll[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))] 
  
  # remove some columns, since FstAll takes too much RAM
  FstAll[,c(1,3:13):= NULL]
  
  # add chi.square test result to each SNP
  if(file.exists("chisq_test_all.csv")){
    chisq_test_all<-fread("chisq_test_all.csv")
    setkey(FstAll,LGn,Pos)
    setkey(chisq_test_all,LGn,Pos)
    FstAll<-FstAll[chisq_test_all]
    rm(chisq_test_all)}


  # Chi squre test for the top 5% pbs.r SNPs: Rob vs. Gos and Rob vs. Say
  if(file.exists("pbs.r_top5.csv")){
    pbs_r_top5<-fread("pbs.r_top5.csv")
  }else{
  chi_fun <- function(v1,v2,v3,v4,v5,v6){
    #v1 is RobN_A, v2 is RobNreads
    r1 <- chisq.test(as.table(rbind(c(v1,v2-v1),c(v3,v4-v3))))$p.value
    r2 <- chisq.test(as.table(rbind(c(v1,v2-v1),c(v5,v6-v5))))$p.value
    return(list(r1,r2))
    }
  
  FstAll[PBS.percentile.unsmoothed > 0.95,c("chi_RG","chi_RS"):=chi_fun(RobN_A,RobNreads,GosN_A,GosNreads,SayN_A,SayNreads),by = seq_len(nrow(FstAll[PBS.percentile.unsmoothed > 0.95]))]
  
  # NA chi square test result indicate no difference AND both lakes in either allele have 0 counts
  FstAll[is.na(chi_RG),chi_RG:=1]
  FstAll[is.na(chi_RS),chi_RS:=1]
  FstAll[,chiOK:= ifelse(chi_RG<0.05 & chi_RS<0.05, T, F)]
  
  pbs_r_top5<-FstAll[PBS.percentile.unsmoothed > 0.95,]
  setorder(pbs_r_top5,-pbs.r)
  fwrite(pbs_r_top5,"pbs.r_top5.csv")
  }
  #To chi-square test all SNPs:
  #FstAll[RobNreads !=0 & GosNreads !=0 & SayNreads !=0,c("chi_RG","chi_RS"):=chi_fun(RobN_A,RobNreads,GosN_A,GosNreads,SayN_A,SayNreads),by = seq_len(nrow(FstAll[RobNreads !=0 & GosNreads !=0 & SayNreads !=0]))]
  #fwrite(FstAll[,.(LG,LGn,Pos,SNP_window,Base_A,Base_a,GosN_A,GosNreads,RobN_A,RobNreads,SayN_A,SayNreads,chi_RG,chi_RS,chiOK)],"chisq_test_all.csv")
  
############ Smoothing ################  
# Smoothed Fst and PBS
rm(list=setdiff(ls(), "FstAll"))

windowsize <- 5000

# Put every SNP into a bin, bin 1000 include all the SNPs from 1-1000
FstAll[,SNP_window:= ((Pos-1)%/%windowsize+1)*windowsize]

# Smooth SNPs, only use the SNPs, which as reads number satisfy conditions
temp_BinSNPs<-FstAll[NreadsOK==T&chiOK==T,.("PBS.nSNPs" = .N,"PBS.r" = mean(pbs.r)),by=.(LGn, SNP_window)]

# Create a data.table with the full length of genome for temp_BinSNPs to join in.
# I should not use chr length from the table "Stickle_chr_lengths.txt", as it truncates a lot of sequenced chromosomes and lost SNPs.
chr <- FstAll[, .(chr_length = max(Pos)), by = LGn]
chr[,cumulative_chrlengths := cumsum(chr_length)]
chr[,cumulative_chrSTART := shift(cumulative_chrlengths,1,0,"lag")]
chr_forjoin<-chr[rep(1:.N,chr_length%/%windowsize+1)][,window_pos := ((1:.N))*windowsize, by = LGn]

# Join the smoothed Fst table with the chr_forjoin, the latter has the full length of chromosome
setkey(chr_forjoin, LGn, window_pos)
setkey(temp_BinSNPs, LGn, SNP_window)

PBS_binned <- temp_BinSNPs[chr_forjoin]
#PBS_binned[,PBS.percentile.score:=1-ecdf(PBS.r)(PBS.r)]
#PBS_binned[,cumulative_window_pos:=SNP_window+cumulative_chrSTART]
setnames(PBS_binned,"SNP_window","Pos")
#PBS_binned_cleaned <- PBS_binned[,!c("chr_length","cumulative_chrlengths","cumulative_chrSTART")]
#BinSNPs[is.na(nSNPs),nSNPs:=0]
#BinSNPs[is.na(PBS.r),PBS.r:=0]

#################################################################
#  4.  Combined genome-wide plots of PoolSeq, Pi, Tajima's D
#################################################################

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

# if bins are not 1000
TajD_r[,binned_pos := windowsize*((Pos-1)%/%windowsize+1),by = LGn ]
TajD_r<-TajD_r[, .(LGn,Pos=binned_pos,D.nSNPs=sum(D.nSNPs,na.rm=T),D.r = mean(D.r,na.rm=T)), by= .(LGn,binned_pos)]

setkey(TajD_r, LGn, Pos)

TajD_binned <- TajD_r[chr_forjoin]
TajD_binned_cleaned <- TajD_binned[,.(LGn,Pos,D.r,D.nSNPs)]

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
setkey(TajD_PBS, LGn, Pos)

# if bins are not greater than 2500
if (windowsize>2500){
  Pi_r[,binned_pos := windowsize*((Pos-1)%/%windowsize+1),by = LGn ]
  Pi_r<-Pi_r[, .(LGn,Pos=binned_pos,Pi.nSNPs=sum(Pi.nSNPs,na.rm=T),Pi.r = mean(Pi.r,na.rm=T)), by= .(LGn,binned_pos)]
  Pi_binned <- Pi_r[chr_forjoin]
}else{
  # rolljoin in data.table can assign values to its closest neighbor. -Inf means backward rolling. Note that Pi.nSNPs is the sum of those 2-3 rows
  Pi_binned <- Pi_r[chr_forjoin,roll = -2500]
}

Pi_binned_cleaned <- Pi_binned[,.(LGn,Pos,Pi.r,Pi.nSNPs)]
setkey(Pi_binned_cleaned,LGn,Pos)


########### Combine all three tables
Pi_TajD_PBS<-TajD_binned_cleaned[PBS_binned]
Pi_TajD_PBS<-Pi_binned_cleaned[Pi_TajD_PBS]



### combine statistics
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
#  5.  Identify main sites of interest
#################################################################

source("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/R Code/stickleback_poolseq_RobGosSay/gwscaR_plot.R")
combined_top5_df <- as.data.frame(combined_top5)
pbs_r_top5_df <- as.data.frame(pbs_r_top5[pbs_r_top5$NreadsOK==T &pbs_r_top5$chiOK==T &PBS.percentile.unsmoothed >0.99,])
PBS_binned_df <- as.data.frame(PBS_binned)
FstAll_df<- as.data.frame(FstAll[NreadsOK==T,c(LGn,Pos,pbs.r)])
Pi_TajD_PBS_cleaned_df <- as.data.frame(Pi_TajD_PBS_cleaned)

fst.plot(fst.dat = PBS_binned_df ,scaffold.widths=chr[,.(LGn,chr_length)],scaffs.to.plot= 21,
        fst.name="pbs.r", chrom.name="LGn", bp.name="Pos",
        y.lim=c(0,6),xlabels=21,xlab.indices=NULL,
        axis.size=0.5,pt.cols=c("darkblue","darkgreen"),pt.cex=0.5)

# Dan's plot function from file Graphics function.R
{
  plotgene <- function(dat, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, popstoplot = c("R", "G", "M")){
  dat[,pos.tr:=Pos/1000000]
  npops <- length(popstoplot)
  par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
  if(specifybps == F){
    minpos <- 0
    maxpos <- dat[LGn==focalLG,max(Pos)]
  }
  ymin <- dat[LGn==focalLG,min(pbs.r,na.rm=T)]
  ymax <- dat[LGn==focalLG,max(pbs.r,na.rm=T)*1.2]
  options(scipen=5)
  if("R" %in% popstoplot) {
    plot(pbs.r~pos.tr,dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,], pch = 16, cex = 0.6, axes = T, xlab = "", ylab = "PBS ROB", col = "darkgreen", ylim = c(ymin,ymax))
    smoothingSpline = smooth.spline(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & !is.na(pbs.r),.(pbs.r,pos.tr)])
    lines(smoothingSpline,col='red')
    #rect(Genestart/1000000, 0,  Geneend/1000000,ymax, border = rgb(1,0,0,0.5) , lwd = 2, col = rgb(1,0,0,0.5))
  }
  #text(Genestart/1000000+0.3, 0.9*ymax, Genename)
  #mtext("Position on LG (mb)", side = 1,  outer = T , line = 1)
  #mtext(paste("LG", focalLG), side = 3,  outer = T )
}

plotgene(PBS_binned,focalLG=1, specifybps = F, minpos, maxpos, Genestart=10000000, Geneend=20000000, Genename="123", popstoplot = c("R"))
}
#################################################################
#  6.  Zoomed-in on sites of interest
#################################################################

gene_loc<-unique(fread("/Users/pengfoen/Documents/Research/Bolnick lab/Analysis_Expression/Data files RAW/GeneLocations.csv"))
gene_name<-fread("/Users/pengfoen/Documents/Research/Bolnick lab/Analysis_Expression/Data files RAW/GeneID_to_GeneName.csv")
gene_LG<-fread("/Users/pengfoen/Documents/Research/Bolnick lab/Analysis_Expression/Data files RAW/GeneID to LinkageGroup.csv",header=T)
setnames(gene_name,c("Ensembl Gene ID","Associated Gene Name"),c("GeneID","gene.name"))
gene_LG[substr(LinkageGroup, 1, 2) == "MT", LGn := 0]
gene_LG[substr(LinkageGroup, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LinkageGroup, 6, nchar(LinkageGroup))))]
gene_LG[substr(LinkageGroup, 1, 2) == "sc", LGn := as.numeric(substr(LinkageGroup, 10, nchar(LinkageGroup)))]
gene_LG[,GeneID:=substring(V1,2)]

setkey(gene_loc,GeneID)
setkey(gene_LG, GeneID)
setkey(gene_name,GeneID)

# lookup all the exons to find the minimum and maximum bp in every gene, and add buffer zone to take into account of the cis-regulatory areas.
gene_loc_full <- gene_loc[, .(start=min(Start)-500,stop=max(Stop)+500), by=GeneID]

# merge all three tables to gene_info
gene_info<-gene_LG[gene_loc_full]
gene_info<-gene_info[gene_name]
gene_info<-gene_info[,!c("V1", "LinkageGroup")]
setkey(gene_info,LGn,start,stop)

# prepare SNPs table
SNP_info<-PBS_binned[PBS.percentile.score<0.05,.(LGn,start=SNP_window-windowsize+1,stop=SNP_window,PBS.r,PBS.nSNPs,PBS.percentile.score)]

SNP_info<-Pi_TajD_PBS[,c("start","stop"):=list(Pos-windowsize+1,Pos)]

setorder(SNP_info,PBS.percentile.score)

SNP_gene_map<-foverlaps(SNP_info,gene_info,by.x = c("LGn","start","stop"))
SNP_gene_map[duplicated(SNP_gene_map,by=c("gene.name"),)][1:20]
fwrite(SNP_gene_map,"SNP_gene_map_5k_bin")


###

SNP_gene_map <- fread("./allele_freq_Fst/snp_gene_map", header = F, sep = "\t")
