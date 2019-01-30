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
  setwd("D:/Foen Peng/OneDrive - University of Connecticut/")
  Pi_TajD_PBS_gene<-fread("./Analyses_Poolseq/Results/Foen/Pi_TajD_PBS_gene.csv")
  focal_filelist = list.files(path = "./Analyses_QTL/result/Foen", pattern="focal_region", full.names = T) 
  focal_QTL <- do.call(rbind,lapply(focal_filelist,function(i)  {cbind(fread(i, header = T), name=strsplit(basename(i),"_focal_region")[[1]][1])})) 
}

#################################################################
#  4.  Load Pi and Tajima's D of Rob population
#################################################################
{
  focal_qtl_combined<-focal_QTL[, as.data.table(reduce(IRanges(start,end),min.gapwidth = 1000000)),by = LGn]
  setnames(focal_qtl_combined,2:3,c("combined_qtl_start","combined_qtl_end"))
  Pi_TajD_PBS_gene<-Pi_TajD_PBS_gene[,c("gene_start","gene_end"):=list(start,stop)]
  focal_qtl_combined_gene<-focal_qtl_combined[Pi_TajD_PBS_gene,.(i.gene_start,i.gene_end),
                                              nomatch = 0L,
                                              on = c("combined_qtl_start<=start","combined_qtl_end<=stop","LGn"),
                                              allow.cartesian = TRUE][]
  test<-Pi_TajD_PBS_gene[]
  
  
}
