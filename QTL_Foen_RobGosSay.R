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
  qtl_all_focal_snp<-qtl_all[qtl_all_focal_region, 
          on = c("Pos.join>start","Pos.join<end","LGn","qtl.trait"),
          nomatch = 0L]
  
  # add the few significant scaffolds to the joined snp data
  qtl_sig_scaffolds<-qtl_all[substr(V1, 1, 2) == "sc" & sig == TRUE,]
  qtl_all_focal_snp<- merge(qtl_sig_scaffolds,qtl_all_focal_snp,all=TRUE)
  
  # some data clearning and saving the file
  setnames(qtl_all_focal_snp,c("Pos.join","Pos.join.1","width"), c("qtl.focal.region.start","qtl.focal.region.end","qtl.focal.region.width"))
  qtl_all_focal_snp[,c("focal.qtl.start","focal.qtl.end"):=NULL]
  fwrite(qtl_all_focal_snp,"./result/Foen/QTL_AllTraits_snp_in_focal_region.csv")
}  

# plot qtl for each trait
{
  for(i in qtl_all_focal_snp[,unqiue(qtl.trait)]){
    qtl_plot<-qtl_all_focal_snp[qtl.trait==i,]
    
    png(sprintf("./result/Foen/%s_QTL.png",i), width = 4000, height = 2000,res=300)
    setkey(qtl_plot,LGn,Pos)
    chrcol <- 2+2*as.numeric(qtl_plot[,LGn] %%2 == 1)
    chrpch <- ifelse(qtl_plot[,sig], 16, 1)
    plot(qtl_plot[,Cumulative_position],qtl_plot[,abs_Z], col=chrcol, type = 'o', axes=F, xlab = "", ylab = "Worm Presence Z_SNP",pch=chrpch, ylim = c(0, 5))
    text(qtl_plot[sig==T,Cumulative_position],qtl_plot[sig==T,abs_Z], labels = qtl_plot[sig==T,snp.id],cex=0.7, pos=2)
    abline(h=qtl_plot[sig==FALSE, max(abs_Z)], lty=5, col='black')
    Chr.mid <- Cumulative_ChrSTART + Chrlengths / 2
    axis(2)
    axis(1, at= Chr.mid, labels = 1:21)
    rect(Cumulative_ChrSTART,0,Cumulative_Chrlengths,5,border = rgb(0,0,0,0.5), lwd = 2)
    rect(qtl_plot_focal[,cum.start], 0,  qtl_plot_focal[,cum.end],5, border = rgb(1,0,0,0.2) , lwd = 2, col = rgb(1,0,0,0.2))
    dev.off()
    
  }
}
