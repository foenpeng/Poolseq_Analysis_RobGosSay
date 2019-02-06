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
  focal_qtl_snp<-fread("./Analyses_QTL/result/Foen/QTL_AllTraits_snp_focal.csv")
  focal_qtl<-na.omit(focal_qtl_snp[,.(qtl.focal.region.end=qtl.focal.region.end[1]),by=c("qtl.trait","LGn","qtl.focal.region.start")])
  #focal_filelist = list.files(path = "./Analyses_QTL/result/Foen", pattern="focal_region", full.names = T) 
  #focal_qtl <- do.call(rbind,lapply(focal_filelist,function(i)  {cbind(fread(i, header = T), name=strsplit(basename(i),"_focal_region")[[1]][1])})) 
}

#################################################################
#  3.  Plot genomic map for each QTL
#################################################################
{
  ##### Generate a table containing all the genes fall to any QTL
  Pi_TajD_PBS_gene[,c("gene_join_start","gene_join_end"):=list(start,stop)]
  focal_qtl[,c("qtl_join_start","qtl_join_end"):=list(qtl.focal.region.start,qtl.focal.region.end)]
  focal_qtl_gene<-focal_qtl[Pi_TajD_PBS_gene,
                            nomatch = 0L,
                            on = c("qtl_join_start<=gene_join_start","qtl_join_end>=gene_join_end","LGn"),
                            allow.cartesian = TRUE]
  focal_qtl_gene[,pbs.focal.gene := FALSE]
  focal_qtl_gene[pbs.r.perc>0.95 & dp.GR.perc>0.95,pbs.focal.gene := TRUE]
  focal_qtl_gene_significant<-focal_qtl_gene[pbs.focal.gene ==TRUE,.SD,by=c("qtl.trait","LGn","GeneID")]
  #fwrite(focal_qtl_gene,"./Analyses_Combined/focal_qtl_gene_significant.csv")
}


{
  plot.qtl <- function(data_gene, data_qtl, data_qtl_snp, statstoplot, stats95,qtl_sig_level){
   
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    gene_buffer<-1e4
    dat<-data_gene[LGn==data_qtl[,LGn] & qtl.trait==data_qtl[,qtl.trait] & qtl.focal.region.start==data_qtl[,qtl.focal.region.start],]
    gene_focal_merge<-dat[pbs.focal.gene==TRUE, as.data.table(reduce(IRanges(start-gene_buffer,stop+gene_buffer),min.gapwidth = 0L)),by = LGn]
    gene_pch <- ifelse(dat[,pbs.focal.gene], 16, 1)
    gene_cex<- ifelse(dat[,pbs.focal.gene], 1.2, 0.5)
    
    qtl_snp_plot<-data_qtl[data_qtl_snp,.(qtl.focal.region.start, qtl.focal.region.end, Pos, abs_Z, snp.id, sig),
                           on=c("qtl.trait","LGn","qtl.focal.region.start"),nomatch=0L]
    #cols_to_divide<-c("qtl.focal.region.start", "qtl.focal.region.end", "Pos")
    #qtl_snp_plot[,(cols_to_divide):=lapply(.SD,function(x) x/(1e6)),.SDcols = cols_to_divide ]
    qtl_snp_pch<-ifelse(qtl_snp_plot[,sig], 16, 1)
    qtl_snp_cex<-ifelse(qtl_snp_plot[,sig], 1.2, 0.5)
    
    for(i in 1:length(statstoplot)){
      var_plot<-statstoplot[i]
      if(var_plot == "qtl"){
        plot(qtl_snp_plot[,.(Pos,abs_Z)],pch=qtl_snp_pch, type = 'o', cex = qtl_snp_cex, axes = T, xlim=c(data_qtl[1,qtl.focal.region.start], data_qtl[1,qtl.focal.region.end]), xlab = NULL, ylab = var_plot, col=col[i])
        text(qtl_snp_plot[,.(Pos, abs_Z)], labels = qtl_snp_plot[,snp.id],cex=1, pos=2)
        abline(h=qtl_sig_level[qtl.trait==data_qtl[,qtl.trait]], lty=5, col='black')
      }else{
        var_95<-stats95[,get(var_plot)]
        print(var_95)
        plot(dat[,.(((start+stop)/2),get(var_plot))], pch=gene_pch, cex = gene_cex, axes = T, xlim=c(data_qtl[1,qtl.focal.region.start], data_qtl[1,qtl.focal.region.end]), xlab = NULL, ylab = var_plot, col=col[i])
        abline(h=var_95, lty=5, col='black')
        if(dat[pbs.focal.gene==TRUE,.N]>0){
          text(dat[pbs.focal.gene==TRUE,.(((start+stop)/2),get(var_plot))], labels = dat[pbs.focal.gene==TRUE,paste(substr(GeneID,14,nchar(GeneID)),gene.name)],cex=1, pos=2)
          rect(gene_focal_merge[,start], 0,  gene_focal_merge[,end],5, border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.4))}
      }
    }
    mtext(sprintf("%s Chromosome: %s", dat[,qtl.trait], dat[,LGn]),side=1,line=1, outer = TRUE)
  }
  
  # create a list for quantile 95% of pbs and significance level of qtl marker.
  popgen_statstoplot<-c("pbs.r","dp.GR")
  
  popgen_stats95<-focal_qtl_gene[,lapply(.SD,quantile,probs=0.95,na.rm=TRUE),.SDcols = popgen_statstoplot]
  qtl_sig_level<-focal_qtl_snp[sig==FALSE,.(sig_Z=max(abs_Z)),by=qtl.trait]
    
  plot.qtl(focal_qtl_gene,focal_qtl[1],focal_qtl_snp, c("qtl",popgen_statstoplot),popgen_stats95,qtl_sig_level)
  
  focal_qtl_plot<-focal_qtl[qtl.trait=="Fibrosis",]
  
  for (j in 1:focal_qtl_plot[,.N]) {
    #png(filename=sprintf("/Users/pengfoen/Documents/Research/Bolnick lab/Analyses_Poolseq/Results/Foen/chr%s.pop.divergence.png", chromosome),width = 3000, height = 2000,res=300)
    #dev.off()
    plot.qtl(focal_qtl_gene,focal_qtl[j], focal_qtl_snp, c("qtl",popgen_statstoplot),popgen_stats95,qtl_sig_level)
  } 
}


focal_qtl_combined<-focal_QTL[, as.data.table(reduce(IRanges(start,end),min.gapwidth = 1000000)),by = LGn]