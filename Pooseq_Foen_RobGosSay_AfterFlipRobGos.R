####  Analysis of PoolSeq data for Weber and Steinel ms, speed up with data.table
#### FP


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

setwd("/Users/pengfoen/OneDrive - University of Connecticut/Analyses_Poolseq")
FstAll<-fread("../Poolseq_8pops/result/RobGosSay_Fst.csv")

### Pre-proceesing Fst data
{
  head(FstAll)
  colnames(FstAll) <- c("LGn", "Pos", "Base_A", "Base_a", "GosN_A","RobN_A","SayN_A", "GosNreads",  "RobNreads",  "SayNreads", "GRFst", "GSFst" , "RSFst")
  FstAll[,GosFreqA := FstAll$GosN_A/FstAll$GosNreads]
  FstAll[,RobFreqA := FstAll$RobN_A/FstAll$RobNreads]
  FstAll[,SayFreqA := FstAll$SayN_A/FstAll$SayNreads]
  
  # remove SNPs that are very similar among these three populations. The original SNP call was based on 8 population comparison.
  FstAll<-FstAll[GosFreqA+RobFreqA+SayFreqA>0.05 & GosFreqA+RobFreqA+SayFreqA<2.95,]
  
  # remove SNPs that have too few or too many reads
  # hist(FstAll[,GosNreads],breaks=100)
  min_cov<- 20
  max_cov<-400
  FstAll<-FstAll[GosNreads >=min_cov & RobNreads >=min_cov & SayNreads >=min_cov & GosNreads <=max_cov & RobNreads <=max_cov & SayNreads <=max_cov,]
} 

#################################################################
#  3.  Calculate PBS and other statistics
#################################################################

{
  # Caculate PBS for Rob branch, Gos branch, and Marine branch
  # From Sequencing of Fifty Human Exomes Reveals Adaptation to High Altitude Supplement
  FstAll[,pbs.r := ((-log(1-FstAll$GRFst )) + (- log(1-FstAll$RSFst )) - (- log(1-FstAll$GSFst )))/2]
  FstAll[,pbs.g := ((-log(1-FstAll$GRFst )) + (- log(1-FstAll$GSFst )) - (- log(1-FstAll$RSFst )))/2]
  FstAll[,pbs.s := ((-log(1-FstAll$RSFst )) + (- log(1-FstAll$GSFst )) - (- log(1-FstAll$GRFst )))/2]
  #FstAll[,pbs.r.percentile:=ecdf(pbs.r)(pbs.r)]
  
  # Calculate dxy between Gos and Rob fishes.dxy = (p1∗(1−p2))+(p2∗(1−p1)). According to "Comparative analysis examining patterns of genomic differentiation across multiple episodes of population divergence in birds" 
  # FstAll[,dxy.GR := GosFreqA  * (1- RobFreqA)  + RobFreqA * (1- GosFreqA)] 
  # FstAll[,dxy.RS := SayFreqA  * (1- RobFreqA)  + RobFreqA * (1- SayFreqA)] 
  # FstAll[,dxy.GS := GosFreqA  * (1- SayFreqA)  + SayFreqA * (1- GosFreqA)] 
  
  # Calculate the absolute difference of allele frequencies
  FstAll[,AFD.GR := abs(GosFreqA  - RobFreqA)] 
  
  # clean some table and columns
  #FstAll[,c(1,3:13,24):= NULL]
  rm(list=setdiff(ls(), "FstAll"))
  
  # save a copy of SNP information
  #FstSave<-FstAll[NreadsOK == T,]
  #fwrite(FstSave,"popgen_info_filter_seq_depth.csv")
  #fwrite(FstAll[NreadsOK==F,],"popgen_info_filter_seq_depth_FALSE.csv")
  
}  

#################################################################
#  4. Binning SNPs into windows for Sliding window analysis
#################################################################

{
  
  chr <- FstAll[, .(chr_snp=.N,chr_length = max(Pos)), keyby = LGn]
  chr[,Cumulative_Chrlengths:=cumsum(chr_length)]
  chr[,Cumulative_ChrSTART := c(0, Cumulative_Chrlengths[-length(Cumulative_Chrlengths)])]
  
  {  windowsize_bp <- 20000
    stepsize_bp<-10000
  
    chr_slide_byBP<-chr[rep(1:.N,(chr_length)%/%stepsize_bp)][,.(window_start=(0:(.N-1))*stepsize_bp, window_end=((0:(.N-1))*stepsize_bp+windowsize_bp),chr_length=chr_length, cum_chr_start=Cumulative_ChrSTART), by = LGn]
    chr_slide_byBP<-chr_slide_byBP[window_end<chr_length+stepsize_bp]
  
    FstAll[,Pos_join:=Pos]
    PBS_slide_bp<-FstAll[chr_slide_byBP, 
                      on=c("LGn","Pos_join>window_start","Pos_join<window_end"), 
                      allow.cartesian=T][,.("PBS.nSNPs" = .N,
                                            "PBS.r" = mean(pbs.r,na.rm=T), 
                                            "PBS.g" = mean(pbs.g,na.rm=T), 
                                            "PBS.s" = mean(pbs.s,na.rm=T), 
                                            "GRFst" = mean(GRFst,na.rm=T), 
                                           # "dxy.GR"= mean(dxy.GR,na.rm=T),
                                           # "dxy.GS"= mean(dxy.GS,na.rm=T),
                                           # "dxy.RS"= mean(dxy.RS,na.rm=T),
                                            "AFD.GR"= mean(AFD.GR,na.rm=T),
                                            "window_end" = Pos_join.1[1],
                                           "Pos_cum" = Pos_join[1] + cum_chr_start[1]),
                                         keyby=.(LGn,Pos_join)]
    setnames(PBS_slide_bp, "Pos_join", "window_start")
    PBS_slide_bp[,Pos:=(window_end+window_start)/2]
    fwrite(PBS_slide_bp, "./Results/Foen/PBS_slide_20kb.csv")
  }

  
  { windowsize_snp <- 500
    stepsize_snp<-250
    
    chr_slide_bySNP<-chr[rep(1:.N,(chr_snp)%/%stepsize_snp)][,.(window_start=(0:(.N-1))*stepsize_snp, window_end=((0:(.N-1))*stepsize_snp+windowsize_snp),chr_snp=chr_snp,chr_length, cum_chr_start=Cumulative_ChrSTART), by = LGn]
    chr_slide_bySNP<-chr_slide_bySNP[window_end<chr_snp+stepsize_snp]
    
    FstAll[,rowID:=rowid(LGn)]
    PBS_slide_snp<-FstAll[chr_slide_bySNP, 
                         on=c("LGn","rowID>=window_start","rowID<window_end"), 
                         allow.cartesian=T][,.("PBS.r" = mean(pbs.r,na.rm=T), 
                                               "PBS.g" = mean(pbs.g,na.rm=T), 
                                               "PBS.s" = mean(pbs.s,na.rm=T), 
                                               "GRFst" = mean(GRFst,na.rm=T), 
                                               # "dxy.GR"= mean(dxy.GR,na.rm=T),
                                               # "dxy.GS"= mean(dxy.GS,na.rm=T),
                                               # "dxy.RS"= mean(dxy.RS,na.rm=T),
                                               "AFD.GR"= mean(AFD.GR,na.rm=T),
                                               "Pos_start" = min(Pos),
                                               "Pos_end" = max(Pos),
                                               "Pos" = (max(Pos)+min(Pos))/2,
                                               "Pos_cum" = (max(Pos)+min(Pos))/2 + max(cum_chr_start)),
                                            keyby=.(LGn,rowID)]
    setnames(PBS_slide_snp, "rowID", "SNP_number")
    fwrite(PBS_slide_snp, "./Results/Foen/PBS_slide_500snp.csv")
    
  }  
}
#################################################################
#  6.  Calculate statistics for every gene 
#################################################################

##### Load gene information
{
  gene_dir <- paste(dirname(getwd()),"/Analysis_Expression/Data files RAW",sep='')
  gene_loc<-unique(fread(paste(gene_dir,"/GeneLocations.csv",sep="")))
  gene_name<-fread(paste(gene_dir,"/GeneID_to_GeneName.csv",sep=""))
  gene_LG<-fread(paste(gene_dir,"/GeneID to LinkageGroup.csv",sep=""),header = TRUE)
  setnames(gene_name,c("Ensembl Gene ID","Associated Gene Name"),c("GeneID","gene.name"))
  gene_LG[substr(LinkageGroup, 1, 2) == "MT", LGn := 0]
  gene_LG[substr(LinkageGroup, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LinkageGroup, 6, nchar(LinkageGroup))))]
  gene_LG[substr(LinkageGroup, 1, 2) == "sc", LGn := as.numeric(substr(LinkageGroup, 10, nchar(LinkageGroup)))]
  gene_LG[,GeneID:=substring(V1,2)]
  
  setkey(gene_loc,GeneID)
  setkey(gene_LG, GeneID)
  setkey(gene_name,GeneID)
  
  # lookup all the exons to find the minimum and maximum position in every gene, and add buffer zone to take into account of the cis-regulatory areas.
  gene_loc_full <- gene_loc[, .(start=min(Start)-500,stop=max(Stop)+500), by=GeneID]
  
  # merge three tables to gene_info
  gene_info<-gene_LG[gene_loc_full]
  gene_info<-gene_info[gene_name]
  gene_info<-gene_info[,!c("V1", "LinkageGroup")]
  setkey(gene_info,LGn,start,stop)
}  

##### calculate statistics for every gene
{
  PBS_gene<-FstAll[gene_info, 
                   on = c("Pos>=start","Pos<=stop","LGn"), 
                   nomatch = NA, 
                   allow.cartesian = TRUE][, .(start = Pos[1],
                                               stop = Pos.1[1],
                                               gene.name = gene.name[1],
                                               gene.length = Pos.1[1]-Pos[1],
                                               snp.count = .N,
                                               pbs.r = mean(pbs.r,na.rm=T),
                                               pbs.g = mean(pbs.g,na.rm=T),
                                               pbs.s = mean(pbs.s,na.rm=T),
                                               AFD.GR = mean(AFD.GR,na.rm=T),
                                               GRFst = mean(GRFst,na.rm=T)),
                                           keyby = .(LGn,GeneID)]
  
  ecdf_fun <- function(x) ecdf(x)(x)
  PBS_gene[,c("pbs.r.perc","pbs.g.perc","pbs.s.perc","AFD.GR.perc","GRFst.perc") := lapply(.SD, ecdf_fun), .SDcols = 8:12]
  
  
  fwrite(PBS_gene, "./Results/Foen/PBS_gene.csv")
}
#################################################################
#  7.  Plotting results
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

##### Plot genomic map for statistics calculated by gene

plotgene_gene <- function(dat, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, statstoplot){
  npops <- length(statstoplot)
  par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
  if(specifybps == F){
    minpos <- 0
    maxpos <- dat[LGn==focalLG,max(stop)+10,000]
  }
  options(scipen=5)
  col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
  for(i in 1:length(statstoplot)){
    
    plot(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,.(((stop-start)/2)/1000000,eval(parse(text = statstoplot[i])))], pch = 16, cex = 0.6, axes = T, xlab = paste("Chromosome",focalLG), ylab = statstoplot[i],col=col[i])
    abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),((stop-start)/2)/1000000],col=rgb(red=0,green=0,blue=0,alpha=100,maxColorValue = 255))
  }
  mtext(paste("Chromosome",focalLG,", every dot represents a gene"),side=1,line=1, outer = TRUE)
}





pbs_r_cutoff <- quantile(Pi_TajD_PBS[,PBS.r],0.95,na.rm=T)
FstAll_gene[,focal:= ifelse(pbs.r> pbs.r.cutoff, T, F)]


for (chromosome in 1:21) {
  png(filename=sprintf("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/Results/Foen/chr%s.pop.divergence.png", chromosome),width = 3000, height = 2000,res=300)
  plotgene(Pi_TajD_PBS,focalLG=chromosome, specifybps = F, minpos, maxpos, statstoplot = c("PBS.r", "AFD.GR", "Ra.r","D.r","Pi.r"))
  dev.off()
} 



##### Plot binned genomic map, modified from Dan's plot function from file Graphics function.R
{ 
  plot_bin_genome <- function(dat, statstoplot){
    par(mfrow = c(length(statstoplot),1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    chrcol <- ifelse(dat[,eval(parse(text = paste(statstoplot,".focal",sep="")))], 1, "grey60")
    for(i in 1:length(statstoplot)){
      y_lim<-c(dat[,min(eval(parse(text = statstoplot[i])),na.rm=T)*1.2],dat[,max(eval(parse(text = statstoplot[i])),na.rm=T)*1.2])
      plot(dat[,Pos_cum],dat[,eval(parse(text = statstoplot[i]))], col=col[i],axes=F, xlab = "Linkage Group", ylab = statstoplot[i],pch=16,cex=0.3, ylim = y_lim)
      #text(dat[sig==T,Cumulative_position],dat[sig==T,abs_Z], labels = dat[sig==T,snp.id],cex=0.7, pos=2)
      abline(h=dat[eval(parse(text = paste(statstoplot[i],".focal",sep="")))==FALSE, max(eval(parse(text = statstoplot[i])))], lty=5, col='black')
      #abline(v=dat[eval(parse(text = paste(statstoplot[i],".focal",sep="")))==TRUE, Pos_cum], lty=1, col=rgb(0,0,0,alpha=0.3))
      Chr.mid <- dat[,(min(Pos_cum)+max(Pos_cum))/2,by=LGn][,V1]
      axis(2)
      axis(1, at= c(dat[,(min(Pos_cum)+max(Pos_cum))/2,keyby=LGn][1:22,V1]), labels = c(1:21,"Un"))
      rect(dat[,min(Pos_cum),keyby=LGn][,V1],y_lim[1],dat[,max(Pos_cum),by=LGn][,V1],y_lim[2],border = NA, lwd = 0.1, col= rep(c(rgb(0,0,0,alpha=0.1),NA),1) )
    }
  }
  
  plot_LG <- function(dat, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, statstoplot){
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    if(specifybps == F){
      minpos <- 0
      maxpos <- dat[LGn==focalLG,max(Pos)]
    }
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    for(i in 1:length(statstoplot)){
      
      plot(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,.(Pos/1000000,eval(parse(text = statstoplot[i])))],pch = 16, cex = 0.6, axes = F, xlab = paste("Chromosome",focalLG), ylab = statstoplot[i],col=col[i])
      abline(h=dat[eval(parse(text = paste(statstoplot[i],".focal",sep="")))==FALSE, max(eval(parse(text = statstoplot[i])))], lty=5, col='black')
      abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),Pos/1000000],col=rgb(red=0,green=0,blue=0,alpha=50,maxColorValue = 255))
      ticks<-seq(0,max(dat[LGn %in% c(focalLG),Pos])/1e6,1)
      axis(1, at=seq(0,maxpos/1e6,0.5), labels = F)
      axis(1, at=seq(0,maxpos/1e6,1), labels = seq(0,maxpos/1e6,1))
      axis(2)
    }
    mtext(paste("Chromosome",focalLG,"(Mb)"),side=1,line=1, outer = TRUE)
  }
  
  cutoff_calculation<-function(var, data=dat_plot, threshold=0.99){
    var.cutoff <- quantile(data[,get(var)],threshold,na.rm=T)
    data[,paste0(var,".focal"):= ifelse(get(var)> var.cutoff, T, F)]
  }
  
  dat_plot <- PBS_slide_snp  
  statstoplot = c("PBS.r", "PBS.g", "PBS.s","GRFst","AFD.GR")
  lapply(statstoplot,cutoff_calculation)
  
  png(sprintf("./Results/Foen/Genome_500snp_window.png",statstoplot), width = 5000, height = 3000,res=600)  
  plot_bin_genome(dat=dat_plot,statstoplot)
  dev.off()

  for (chromosome in 1:21) {
    png(filename=sprintf("./Results/Foen/chr_map_20kb_slide/chr%s.pop.divergence_0.99.png", chromosome),width = 3000, height = 2000,res=300)
    plot_LG(dat_plot,focalLG=chromosome, specifybps = F, minpos, maxpos, statstoplot = c("PBS.r", "PBS.g", "PBS.s","GRFst","AFD.GR"))
    dev.off()
  } 
  
}
#################################################################
#  8.  Zoomed-in on sites of interest
#################################################################



# prepare SNPs table
SNP_info<-PBS_binned[PBS.percentile.score<0.05,.(LGn,start=SNP_window-windowsize+1,stop=SNP_window,PBS.r,PBS.nSNPs,PBS.percentile.score)]

SNP_info<-Pi_TajD_PBS[,c("start","stop"):=list(Pos-windowsize+1,Pos)]

setorder(SNP_info,PBS.percentile.score)

SNP_gene_map<-foverlaps(SNP_info,gene_info,by.x = c("LGn","start","stop"))
SNP_gene_map[duplicated(SNP_gene_map,by=c("gene.name"),)][1:20]
fwrite(SNP_gene_map,"SNP_gene_map_5k_bin")


###

SNP_gene_map <- fread("./allele_freq_Fst/snp_gene_map", header = F, sep = "\t")
