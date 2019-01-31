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
  setwd("/Users/pengfoen/OneDrive - University of Connecticut/")
  Pi_TajD_PBS_gene<-fread("./Analyses_Poolseq/Results/Foen/Pi_TajD_PBS_gene.csv")
  focal_filelist = list.files(path = "./Analyses_QTL/result/Foen", pattern="focal_region", full.names = T) 
  focal_qtl <- do.call(rbind,lapply(focal_filelist,function(i)  {cbind(fread(i, header = T), name=strsplit(basename(i),"_focal_region")[[1]][1])})) 
}

#################################################################
#  3.  Plot genomic map for each QTL
#################################################################
{
  ##### Generate a table containing all the genes fall to any QTL
  setnames(focal_qtl,2:3,c("qtl_start","qtl_end"))
  Pi_TajD_PBS_gene[,c("gene_join_start","gene_join_end"):=list(start,stop)]
  focal_qtl[,c("qtl_join_start","qtl_join_end"):=list(qtl_start,qtl_end)]
  focal_qtl_gene<-focal_qtl[Pi_TajD_PBS_gene,
                            nomatch = 0L,
                            on = c("qtl_join_start<=gene_join_start","qtl_join_end>=gene_join_end","LGn"),
                            allow.cartesian = TRUE]
  #fwrite(focal_qtl_gene,"./Analyses_Combined/focal_qtl_gene.csv")
}


{
  plot.qtl <- function(data_gene, data_qtl, statstoplot, stats95){
    print(stats95)
   
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    options(scipen=5)
    gene_buffer<-1e4
    dat<-data_gene[LGn==data_qtl[,LGn] & name==data_qtl[,name] & qtl_start==data_qtl[,qtl_start],]
    gene_focal_merge<-dat[pbs.focal.gene==TRUE, as.data.table(reduce(IRanges(start-gene_buffer,stop+gene_buffer),min.gapwidth = 0L)),by = LGn]
    gene_pch <- ifelse(dat[,pbs.focal.gene], 16, 1)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    for(i in 1:length(statstoplot)){
      var_plot<-statstoplot[i]
      var_95<-stats95[i]
      if(var_plot != "qtl"){
        plot(dat[,.(((start+stop)/2),get(var_plot))], pch=gene_pch, cex = 0.6, axes = T, xlim=c(data_qtl[1,qtl_start], data_qtl[1,qtl_end]), xlab = NULL, ylab = var_plot, col=col[i])
        abline(h=var_95, lty=5, col='black')
        if(dat[pbs.focal.gene==TRUE,.N]>0){
          text(dat[pbs.focal.gene==TRUE,.(((start+stop)/2),get(var_plot))], labels = dat[pbs.focal.gene==TRUE,paste(substr(GeneID,14,nchar(GeneID)),gene.name)],cex=0.7, pos=2)
          rect(gene_focal_merge[,start], 0,  gene_focal_merge[,end],5, border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.4))}
      }
    }
    mtext(sprintf("%s Chromosome: %s", dat[,name], dat[,LGn]),side=1,line=1, outer = TRUE)
  }
  
  # create a list for quantile 95%
  focal_qtl_gene[,pbs.focal.gene := FALSE]
  focal_qtl_gene[pbs.r.perc>0.95 & dp.GR.perc>0.95,pbs.focal.gene := TRUE]
  statstoplot=c("pbs.r","dp.GR")
  stats95=do.call(rbind,lapply(focal_qtl_gene[,statstoplot,with=FALSE],quantile,probs=0.95,na.rm=TRUE))
  plot.qtl(focal_qtl_gene,focal_qtl[1], statstoplot,stats95)
  focal_qtl_plot<-focal_qtl[name=="QTL_fibrosis"& LGn==2,]
  for (j in 1:focal_qtl_plot[,.N]) {
    #png(filename=sprintf("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/Results/Foen/chr%s.pop.divergence.png", chromosome),width = 3000, height = 2000,res=300)
    #dev.off()
    plot.qtl(focal_qtl_gene,focal_qtl[j], statstoplot,stats95)
  } 
}


focal_qtl_combined<-focal_QTL[, as.data.table(reduce(IRanges(start,end),min.gapwidth = 1000000)),by = LGn]