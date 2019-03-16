#################################################################
#  1.  Load packages
#################################################################
{
  library(data.table)
  library(IRanges)
}

#################################################################
#  2.  Load data
#################################################################
{
  rm(list=ls())
  setwd("/Users/pengfoen/OneDrive - University of Connecticut/Analyses_QTL")
  setwd("D:/Foen Peng/OneDrive - University of Connecticut/Analyses_QTL")
  qtl_singlemarker_filelist = list.files(path = "./result/SingleMarkerResultsDan/", pattern="_SingleMarkerResults.csv", full.names = T) 
  qtl_singlemarker <- do.call(rbind,lapply(qtl_singlemarker_filelist[-5],function(i)  {cbind(fread(i, header = T)[,2:7], qtl.trait=strsplit(basename(i),"_SingleMarkerResults.csv")[[1]][1])})) 
  
  qtl_snp_pos<- fread("./QTL_snp_pos_corrected.csv", header=TRUE)
  qtl_all<-qtl_snp_pos[qtl_singlemarker, on = "snpnum"]
  rm(list=setdiff(ls(), "qtl_all"))
  
  #qtl_snp_pos[substr(V1, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(V1, 6, nchar(V1))))]
  #qtl_snp_pos[substr(V1, 1, 2) == "sc", LGn := as.numeric(substr(V1, 10, nchar(V1)))]
  #snp_pos_correct<-fread("./result/snp_pos_correct.csv", header=TRUE)
  #snp_pos_correct1<-qtl_snp_pos[snp_pos_correct, 
   #                             on = c("LGn==lg","V2==snp.id"), 
    #                            nomatch = NA]
  #fwrite(snp_pos_correct1, "QTL_snp_pos_corrected.csv")
}
#################################################################
#  3.  Plot QTL map
#################################################################

# Define a continuous genomic x axis
{
  Chrlengths <- read.table("Stickle_chr_lengths.txt", sep = ",")
  Chrlengths <- as.vector(as.matrix(Chrlengths) )
  Cumulative_Chrlengths <- cumsum(Chrlengths)
  Cumulative_ChrSTART <- c(0, Cumulative_Chrlengths[-length(Cumulative_Chrlengths)])
  chr_map_full<-setDT(list(LGn=1:21,
                           Chrlengths=Chrlengths,
                           Cumulative_Chrlengths=Cumulative_Chrlengths,
                           Cumulative_ChrSTART=Cumulative_ChrSTART))
}

# Process the qtl data
{
  # add 2 Mb buffer to each end of all significant markers
  qtl_buffer<-2000000
  qtl_all[,sig:=(P_SNP<0.05)]
  qtl_all[,abs_Z:=abs(Z_SNP)]
  # I put the snp marker with na position in the beginning of the LGn.
  qtl_all[,Cumulative_position:=ifelse(is.na(Pos), Cumulative_ChrSTART[LGn], Pos + Cumulative_ChrSTART[LGn])]
  qtl_all[sig==TRUE, c("focal.qtl.start", "focal.qtl.end"):=list(ifelse(Pos-qtl_buffer<=0,0,Pos-qtl_buffer),
                                                                 ifelse(Pos+qtl_buffer>=Chrlengths[LGn],Chrlengths[LGn],Pos+qtl_buffer))]
  # Merge the focal qtl regions which have some overlaps.
  qtl_all_focal_region<-qtl_all[sig==TRUE & is.na(Pos)==FALSE, as.data.table(reduce(IRanges(focal.qtl.start,focal.qtl.end),min.gapwidth = 0L)),by = .(LGn,qtl.trait)]
  qtl_all_focal_region[,c("qtl.focal.region.cum.start","qtl.focal.region.cum.end"):=list(start+Cumulative_ChrSTART[LGn],end+Cumulative_ChrSTART[LGn])]
  
  # join the focal region with the snp data
  qtl_all[,Pos.join:=Pos]
  qtl_all_focal_region[,c("qtl.focal.region.start","qtl.focal.region.end"):=.(start,end)]
  qtl_all_focal_snp<-qtl_all_focal_region[qtl_all, 
          on = c("start<Pos.join","end>Pos.join","LGn","qtl.trait"),
          nomatch = NA]

  # some data clearning and saving the file
  setnames(qtl_all_focal_snp,"width", "qtl.focal.region.width")
  qtl_all_focal_snp[,c("focal.qtl.start","focal.qtl.end","start","end"):=NULL]
  fwrite(qtl_all_focal_snp,"./result/Foen/QTL_AllTraits_snp_in_focal_region.csv")
}  

# plot qtl for each trait
{
  for(i in qtl_all_focal_snp[,unique(qtl.trait)]){
    qtl_plot<-subset(qtl_all_focal_snp, qtl.trait==i)
    
    png(sprintf("./result/Foen/%s_QTL.png",i), width = 4000, height = 2000,res=300)
    setkey(qtl_plot,LGn,Pos)
    chrcol <- 2+2*as.numeric(qtl_plot[,LGn] %%2 == 1)
    chrpch <- ifelse(qtl_plot[,sig], 16, 1)
    plot(qtl_plot[,Cumulative_position],qtl_plot[,abs_Z], col=chrcol, type = 'o', axes=F, xlab = "Linkage Group", ylab = sprintf("%s Z_SNP",i),pch=chrpch, ylim = c(0, 6))
    text(qtl_plot[sig==T,Cumulative_position],qtl_plot[sig==T,abs_Z], labels = qtl_plot[sig==T,snp.id],cex=0.7, pos=2)
    abline(h=qtl_plot[sig==FALSE, max(abs_Z)], lty=5, col='black')
    Chr.mid <- Cumulative_ChrSTART + Chrlengths / 2
    axis(2)
    axis(1, at= Chr.mid, labels = 1:21)
    rect(Cumulative_ChrSTART,0,Cumulative_Chrlengths,6,border = rgb(0,0,0,0.5), lwd = 2)
    focal_region_plot<-na.omit(qtl_plot[,.(unique(qtl.focal.region.cum.start),unique(qtl.focal.region.cum.end))])
    rect(focal_region_plot[,V1], 0,  focal_region_plot[,V2],6, border = rgb(1,0,0,0.2) , lwd = 2, col = rgb(1,0,0,0.2))
    dev.off()
    
  }
}

# plot qtl for each LG

plot_qtl_LGn <- function(dat, LG, maxpos, sig_level,statstoplot){
  setkey(dat,LGn,Pos)
  npops <- length(statstoplot)
  par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
  options(scipen=5)
  col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta","gold4")
  for(i in 1:length(statstoplot)){
    dat_subplot<-dat[qtl.trait==statstoplot[i]]
    chrpch <- ifelse(dat_subplot[,sig], 16, 1)
    plot(dat_subplot[,.(Pos/1e6,abs_Z)],pch = chrpch, cex = 1.5, type="o", ylim=c(0,5), xlim=c(0,maxpos/1e6), axes = F, xlab = paste("Chromosome",LG), ylab = statstoplot[i],col=col[i])
    abline(h=sig_level[qtl.trait==statstoplot[i],V1], lty=5, col='black')
    text(dat_subplot[,.(Pos/1e6,abs_Z)], labels = dat_subplot[,snp.id],cex=1, pos=3)
    if(dat_subplot[sig==T,.N]>0){
      abline(v=dat_subplot[sig==T,Pos/1e6],lwd=4, col=rgb(red=0,green=0,blue=0,alpha=100,maxColorValue = 255))
      text(dat_subplot[sig==T,.(Pos/1e6,abs_Z-2)], labels = dat_subplot[sig==T,round(Pos/1e6,2)],cex=1, pos=3)
    }
    axis(1, at=seq(0,maxpos/1e6,0.5), labels = F)
    axis(1, at=seq(0,maxpos/1e6,1), labels = seq(0,maxpos/1e6,1))
    axis(2)
  }
  mtext(paste("Chromosome",LG),side=1,line=1, outer = TRUE)
}


trait<-unique(qtl_all[,qtl.trait])
sig_level<-qtl_all[sig==FALSE, max(abs_Z), by=qtl.trait]
for(LG in 1:21){
  png(sprintf("./result/Foen/QTL_per_LG/LG%s_QTL.png",LG), width = 3000, height = 2000,res=300)
  plot_qtl_LGn(qtl_all[LGn==LG,],LG,maxpos=Chrlengths[LG], sig_level,statstoplot = trait)
  dev.off()
}
