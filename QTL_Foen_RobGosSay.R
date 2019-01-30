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
  #setwd("/Users/pengfoen/OneDrive - University of Connecticut/Analyses_QTL")
  setwd("D:/Foen Peng/OneDrive - University of Connecticut/Analyses_QTL")
  qtl_fib<-fread("./result/SingleMarkerResultsDan/WormPresent_SingleMarkerResults.csv", header=TRUE)
  qtl_snp_pos<- fread("./QTL_snp_pos_corrected.csv", header=TRUE)
  qtl_fib<-qtl_fib[qtl_snp_pos, on = "snpnum"]
  rm(list=setdiff(ls(), "qtl_fib"))
  
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

##### Plot whole genomic map

# Define a continuous genomic x axis, MTDNA near zero
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

# Plot QTL map 
{
  # pre-processing for plot
  qtl_buffer<-2000000
  qtl_fib[,sig:=(P_SNP<0.05)]
  qtl_fib[,abs_Z:=abs(Z_SNP)]
  qtl_fib[,Cumulative_position:=ifelse(is.na(Pos), Cumulative_ChrSTART[LGn], Pos + Cumulative_ChrSTART[LGn])]
  qtl_fib[sig==TRUE, c("focal.qtl.start", "focal.qtl.end"):=list(ifelse(Pos-qtl_buffer<=0,0,Pos-qtl_buffer),
                                                                 ifelse(Pos+qtl_buffer>=Chrlengths[LGn],Chrlengths[LGn],Pos+qtl_buffer))]
  qtl_fib_focal<-qtl_fib[sig==TRUE & is.na(Pos)==FALSE, as.data.table(reduce(IRanges(focal.qtl.start,focal.qtl.end),min.gapwidth = 0L)),by = LGn]
  qtl_fib_focal[,c("cum.start","cum.end"):=list(start+Cumulative_ChrSTART[LGn],end+Cumulative_ChrSTART[LGn])]
  fwrite(qtl_fib_focal,"./result/Foen/QTL_WormPresence_focal_region.csv")
  
  # plot qtl 
  png("./result/Foen/WormPresence_QTL.png", width = 4000, height = 2000,res=300)
  setkey(qtl_fib,LGn,Pos)
  chrcol <- 2+2*as.numeric(qtl_fib[,LGn] %%2 == 1)
  chrpch <- ifelse(qtl_fib[,sig], 16, 1)
  plot(qtl_fib[,Cumulative_position],qtl_fib[,abs_Z], col=chrcol, type = 'o', axes=F, xlab = "", ylab = "Worm Presence Z_SNP",pch=chrpch, ylim = c(0, 5))
  text(qtl_fib[sig==T,Cumulative_position],qtl_fib[sig==T,abs_Z], labels = qtl_fib[sig==T,snp.id],cex=0.7, pos=2)
  abline(h=qtl_fib[sig==FALSE, max(abs_Z)], lty=5, col='black')
  Chr.mid <- Cumulative_ChrSTART + Chrlengths / 2
  axis(2)
  axis(1, at= Chr.mid, labels = 1:21)
  rect(Cumulative_ChrSTART,0,Cumulative_Chrlengths,5,border = rgb(0,0,0,0.5), lwd = 2)
  rect(qtl_fib_focal[,cum.start], 0,  qtl_fib_focal[,cum.end],5, border = rgb(1,0,0,0.2) , lwd = 2, col = rgb(1,0,0,0.2))
  dev.off()
}
