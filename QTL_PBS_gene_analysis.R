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
rm(list=ls())

## function to convert linkage group names from string to numeric number
  convert_LG_name<-function(dataset, LG_column){
    dataset[substr(get(LG_column), 1, 2) == "MT", LGn := 0]
    dataset[substr(get(LG_column), 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
    dataset[substr(get(LG_column), 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  }

{

  setwd("D:/Foen Peng/OneDrive - University of Connecticut/")
  setwd("/Users/pengfoen/OneDrive - University of Connecticut/")
  Pi_TajD_PBS_gene<-fread("./Analyses_Poolseq/Results/Foen/Pi_TajD_PBS_gene.csv")
  focal_qtl_snp<-fread("./Analyses_QTL/result/Foen/QTL_AllTraits_snp_focal.csv")
  focal_qtl<-na.omit(focal_qtl_snp[,.(qtl.focal.region.end=qtl.focal.region.end[1]),by=c("qtl.trait","LGn","qtl.focal.region.start")])
  qtl_traits <- sort(unique(focal_qtl$qtl.trait))
  manual_qtl_check<-fread("./Analyses_QTL/manual_qtl_check.csv")
  
  # inport the expression data
  expr_pop<-fread("./Analysis_Expression/Tables of results/Front Immunol LM Results Population.csv",select=1:8)
  expr_pop[,source:="pop"]
  expr_int<-fread("./Analysis_Expression/Tables of results/Front Immunol LM Results Interaction.csv",select=1:8)
  expr_int[,source:="int"]
  expr_infection<-fread("./Analysis_Expression/Tables of results/Front Immunol LM Results Infection.csv",select=1:8)
  expr_infection[,source:="infect"]
  expr<-merge(expr_pop,expr_int, all=TRUE)
  expr<-merge(expr,expr_infection,all=TRUE)
  
  # import the GWAS data
  GWAS_filelist = list.files(path = "./Analyses_GWAS/Genome_scans_summary/tabs", full.names = T) 
  GWAS <- do.call(rbind,lapply(GWAS_filelist,function(i)  {cbind(fread(i,header=T, col.names = c("LG","Pos","Start","End")), name=strsplit(basename(i),"Spreadsheet")[[1]][1])})) 
  GWAS <- na.omit(GWAS, cols="Pos")
  GWAS <- convert_LG_name(GWAS, "LG")
  
  # read popgen files, which is very large
  if(file.exists("./Analyses_Poolseq/popgen_info_filter_seq_depth_chi_pbs.csv")){
    popgen_filter<-fread("./Analyses_Poolseq/popgen_info_filter_seq_depth_chi_pbs.csv",header = T)
  }else{
    popgen_all<-fread("./Analyses_Poolseq/popgen_info_filter_seq_depth.csv",header = T)
    # filter at 95% quantile
    popgen_filter<-popgen_all[chiOK==TRUE & pbs.r>1.4358 & dp.GR > 0.828,]
    fwrite(popgen_filter,"./Analyses_Poolseq/popgen_info_filter_seq_depth_chi_pbs.csv")
    rm(popgen_all)
  }

  #plot(popgen_filt[,.(dp.GR,pbs.r)])
  
  #focal_filelist = list.files(path = "./Analyses_QTL/result/Foen", pattern="focal_region", full.names = T) 
  #focal_qtl <- do.call(rbind,lapply(focal_filelist,function(i)  {cbind(fread(i, header = T), name=strsplit(basename(i),"_focal_region")[[1]][1])})) 
}

#################################################################
#  3.  Prepare gene level average pbs data for genomic map
#################################################################

##### Generate a table containing all the genes fall to any QTL
{
  Pi_TajD_PBS_gene[,c("gene_join_start","gene_join_end"):=list(start,stop)]
  focal_qtl[,c("qtl_join_start","qtl_join_end"):=list(qtl.focal.region.start,qtl.focal.region.end)]
  focal_qtl_gene<-focal_qtl[Pi_TajD_PBS_gene,
                            nomatch = 0L,
                            on = c("qtl_join_start<=gene_join_start","qtl_join_end>=gene_join_end","LGn"),
                            allow.cartesian = TRUE]
  focal_qtl_gene[,avg.sig.gene := FALSE]
  focal_qtl_gene[(pbs.r.perc>0.95 | pbs.g.perc>0.95) & dp.GR.perc>0.95,avg.sig.gene := TRUE]
  }  

##### Save a copy of significant genes without duplicates
{
  focal_qtl_gene_significant<-focal_qtl_gene[avg.sig.gene ==TRUE,.SD,by=c("qtl.trait","LGn","GeneID")]
  focal_qtl_gene_significant[,paste0("In.", qtl_traits):=lapply(qtl_traits, function(i) as.integer(qtl.trait==i))]
  focal_qtl_gene_significant[,(26:31):=lapply(.SD, sum),.SDcols=26:31, by=GeneID]
  focal_qtl_gene_significant<-unique(focal_qtl_gene_significant, by="GeneID")
  focal_qtl_gene_significant[,c("qtl.trait","qtl.focal.region.start","qtl.focal.region.end","qtl_join_start","qtl_join_end","avg.sig.gene"):=NULL]
  #fwrite(focal_qtl_gene_significant,"./Analyses_Combined/results/Foen/focal_qtl_gene_significant.csv")
}

#################################################################
#  4.  Prepare individual snp pbs data for genomic map
#################################################################

##### Generate a table that contains SNP after thresholding and falls in QTL
{
  popgen_filter[,Pos.join:=Pos]
  focal_popgen_filter<-focal_qtl[popgen_filter,
                             nomatch=0L,
                             on = c("qtl_join_start<=Pos.join","qtl_join_end>=Pos.join","LGn"),
                             allow.cartesian=TRUE]
  # join the pogen info into gene info
  setnames(focal_qtl_gene,c("qtl_join_start","qtl_join_end"),c("gene_join_start","gene_join_end"))
  focal_popgen_filter[,Pos.join:=Pos]
  focal_popgen_filter_plot<-focal_qtl_gene[focal_popgen_filter,
                                       nomatch=NA,
                                       on = c("gene_join_start<=Pos.join","gene_join_end>=Pos.join","LGn","qtl.trait"),
                                       allow.cartesian=TRUE][,.(qtl.trait,i.qtl.focal.region.start,i.qtl.focal.region.end,
                                                                GeneID,gene.name,start,stop,gene.length,
                                                                i.pbs.r,i.pbs.g,i.dp.GR),by=.(LGn,Pos)]
  setnames(focal_popgen_filter_plot,c("i.qtl.focal.region.start","i.qtl.focal.region.end"),c("qtl.focal.region.start","qtl.focal.region.end"))
  focal_popgen_filter_plot[,snp.sig:=FALSE]
  focal_popgen_filter_plot[i.dp.GR>0.95,snp.sig:=TRUE]
}

##### Generat a table of the highest value SNPs mapping to gene
{
  focal_popgen_filter_save<-focal_popgen_filter
  focal_popgen_filter_save[,paste0("In.", qtl_traits):=lapply(qtl_traits, function(i) as.integer(qtl.trait==i))]
  
  focal_popgen_filter_save[,(33:38):=lapply(.SD, sum),.SDcols=33:38, by=c("LGn","Pos")]
  focal_popgen_filter_save<-unique(focal_popgen_filter_save, by=c("LGn","Pos"))
  focal_popgen_filter_save[,c("qtl.trait","qtl.focal.region.start","qtl.focal.region.end","qtl_join_start","qtl_join_end","chiOK","NreadsOK"):=NULL]
  focal_popgen_filter_save[,Pos.join:=Pos]
  focal_popgen_Snp_map_to_gene<-Pi_TajD_PBS_gene[focal_popgen_filter_save,
                   on = c("gene_join_start<Pos.join","gene_join_end>=Pos.join","LGn"),
                   nomatch = 0L,
                   allow.cartesian=TRUE][,.(start,stop,gene.name,gene.length,snp.count,pbs.r.perc,pbs.g.perc,dp.GR.perc,
                                            GosFreqA,RobFreqA,SayFreqA,i.Base_A,i.GosN_A,i.GosNreads,i.RobN_A,i.RobNreads,i.SayN_A,i.SayNreads,i.pbs.r,i.pbs.g,i.dp.GR,
                                            In.Fibrosis,In.FibrosisORGranuloma,In.Granuloma,In.MaxWormMass,In.ROSresid,In.WormPresent),by=.(LGn,GeneID,Pos)]
  focal_popgen_Snp_map_to_gene_dpGR95<-focal_popgen_Snp_map_to_gene[i.dp.GR>0.95]
  #fwrite(focal_popgen_Snp_map_to_gene,"./Analyses_Combined/results/Foen/focal_popgen_Snp_map_to_gene.csv")
  #fwrite(focal_popgen_Snp_map_to_gene[i.dp.GR>0.95 & i.pbs.g<1],"./Analyses_Combined/results/Foen/focal_popgen_Snp_map_to_gene_dpGR95.csv")
}

##### Look for overlaps between gene level average candidates and snp level candidates
{
  focal_popgen_Snp_map_to_gene_dpGR95_grouped<-focal_popgen_Snp_map_to_gene_dpGR95[,.(start=start[1],stop=stop[1],gene.name=gene.name[1],gene.length=gene.length[1],snp.count=snp.count[1],pbs.r.perc=pbs.r.perc[1],pbs.g.perc=pbs.g.perc[1],dp.GR.perc=dp.GR.perc[1],
                                                                                      In.Fibrosis=In.Fibrosis[1],In.FibrosisORGranuloma=In.FibrosisORGranuloma[1],In.Granuloma=In.Granuloma[1],In.MaxWormMass=In.MaxWormMass[1],In.ROSresid=In.ROSresid[1],In.WormPresent=In.WormPresent[1],
                                                                                      sig.snp.count=.N), by=.(LGn,GeneID)]
  #fwrite(focal_popgen_Snp_map_to_gene_dpGR95_grouped,"./Analyses_Combined/results/Foen/focal_popgen_Snp_map_to_gene_dpGR95_grouped.csv")
  
  overlap_twolist<-focal_popgen_Snp_map_to_gene_dpGR95_grouped[focal_qtl_gene_significant,on="GeneID",nomatch=0L]
  #fwrite(overlap_twolist,"./Analyses_Combined/results/Foen/overlap_between_avg_and_ind_genelist.csv" )
  all_sig_gene<-merge(focal_popgen_Snp_map_to_gene_dpGR95_grouped,focal_qtl_gene_significant,all=TRUE, by=intersect(colnames(focal_popgen_Snp_map_to_gene_dpGR95_grouped),colnames(focal_qtl_gene_significant)))
}

#################################################################
#  5.  Combine QTL and Popgen results with GWAS and Expression
#################################################################

##### combine with expression
{
  qtl_popgen_expr<-expr[all_sig_gene,on="geneID==GeneID",nomatch=0L]
  #fwrite(qtl_popgen_expr,"./Analyses_Combined/results/Foen/qtl_popgen_expr_overlap.csv")
  
} 
  
##### combine with GWAS
{
  setkey(GWAS, LGn, Start, End)
  qtl_popgen_gwas<-foverlaps(all_sig_gene,GWAS, by.x = c("LGn","start","stop"), by.y = c("LGn","Start","End"),type="any",nomatch=0L)
  qtl_popgen_gwas[,c("Start","End"):=NULL]
  #fwrite(qtl_popgen_gwas,"./Analyses_Combined/results/Foen/qtl_popgen_gwas_overlap.csv")
}

##### genes fall into manual QTL check
  qtl_expr<-manual_qtl_check[qtl_popgen_expr, on = c("manual_check_region_start<start","manual_check_region_end>stop","LGn"),nomatch=0L]
  qtl_gwas<-manual_qtl_check[qtl_popgen_gwas, on = c("manual_check_region_start<start","manual_check_region_end>stop","LGn"),nomatch=0L]
  qtl_popgen<-manual_qtl_check[all_sig_gene,, on = c("manual_check_region_start<start","manual_check_region_end>stop","LGn"),nomatch=0L]
  exist_list<-fread("./Analyses_Combined/results/Foen/candidate_list_all_evidence.csv")
  checked_genes<-qtl_popgen[exist_list,on="GeneID",nomatch=0L]
  checked_geneID<-checked_genes[,GeneID]
  unchecked_genes<-qtl_popgen[!(GeneID %in% checked_geneID),]
  fwrite(unchecked_genes,"./Analyses_Combined/results/Foen/unchecked_genes.csv")
  
#### just check QTL with GWAS and expression
{
  qtl_expr_avg<-focal_qtl_gene[expr,on="GeneID==geneID",nomatch=0L]
  qtl_expr_snp<-focal_popgen_Snp_map_to_gene[expr,on="GeneID==geneID",nomatch=0L]
  qtl_expr_avg<-unique(qtl_expr_avg,by="GeneID")
  qtl_expr_snp<-unique(qtl_expr_snp,by="GeneID")
  qtl_expr<-merge(qtl_expr_avg,
                  qtl_expr_snp[,.(start=start[1],stop=stop[1],gene.name=gene.name[1],gene.length=gene.length[1],snp.count=snp.count[1],pbs.r.perc=pbs.r.perc[1],pbs.g.perc=pbs.g.perc[1],dp.GR.perc=dp.GR.perc[1],
                                               sig.snp.count=.N), by=.(LGn,GeneID)], 
                  by=c("LGn","GeneID"),all=TRUE)
  qtl_expr<-manual_qtl_check[qtl_expr, on = c("manual_check_region_start<gene_join_start","manual_check_region_end>gene_join_end","LGn"),nomatch=0L]
  
}  
  {
  setkey(GWAS, LGn, Start, End)
  GWAS_manualQTL<-foverlaps(manual_qtl_check, GWAS, by.x = c("LGn","manual_check_region_start","manual_check_region_end"), by.y = c("LGn","Start","End"),type="any",nomatch=0L)
  setkey(GWAS_manualQTL, LGn, Start, End)
  GWAS_manualQTL_gene<-foverlaps(focal_qtl_gene,GWAS_manualQTL, by.x = c("LGn","start","stop"), by.y = c("LGn","Start","End"),type="any",nomatch=0L)
}
  
#################################################################
#  6.  Plot genomic map for each QTL
#################################################################

##### Plot function
{
  plot.qtl <- function(data_gene, data_qtl, data_qtl_snp,data_popgen, statstoplot, stats95,qtl_sig_level){
    print(data_qtl)
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta","gold4")
    gene_buffer<-1e4
    dat<-data_gene[LGn==data_qtl[,LGn] & qtl.trait==data_qtl[,qtl.trait] & qtl.focal.region.start==data_qtl[,qtl.focal.region.start],]
    avg.sig.gene.merge<-dat[avg.sig.gene==TRUE, as.data.table(reduce(IRanges(start-gene_buffer,stop+gene_buffer),min.gapwidth = 0L)),by = LGn]
    gene_pch <- ifelse(dat[,avg.sig.gene], 16, 1)
    gene_cex<- ifelse(dat[,avg.sig.gene], 1.2, 0.5)
    
    qtl_snp_plot<-data_qtl[data_qtl_snp,.(qtl.focal.region.start, qtl.focal.region.end, Pos, abs_Z, snp.id, sig),
                           on=c("qtl.trait","LGn","qtl.focal.region.start"),nomatch=0L]
    qtl_snp_pch<-ifelse(qtl_snp_plot[,sig], 16, 1)
    qtl_snp_cex<-ifelse(qtl_snp_plot[,sig], 1.2, 0.5)
    setkey(qtl_snp_plot,Pos)
    
    popgen_snp_plot<-data_qtl[data_popgen,on=c("qtl.trait","LGn","qtl.focal.region.start"),nomatch=0L]
    popgen_label<-popgen_snp_plot[snp.sig==TRUE & is.na(GeneID)==FALSE,.SD[which.max(i.pbs.r)], by=GeneID]
    snp.sig.gene.merge<-popgen_snp_plot[snp.sig==TRUE & is.na(GeneID)==FALSE, as.data.table(reduce(IRanges(start-gene_buffer,stop+gene_buffer),min.gapwidth = 0L)),by = LGn]
    popgen_snp_pch<-ifelse(popgen_snp_plot[,snp.sig],16,1)
    
    for(i in 1:length(statstoplot)){
      var_plot<-statstoplot[i]
      if(var_plot == "qtl"){
        plot(qtl_snp_plot[,.(Pos,abs_Z)],pch=qtl_snp_pch, type = 'o', cex = qtl_snp_cex, axes = F, xlim=c(data_qtl[1,qtl.focal.region.start], data_qtl[1,qtl.focal.region.end]), xlab = NULL, ylab = var_plot, col=col[i])
        text(qtl_snp_plot[,.(Pos, abs_Z)], labels = qtl_snp_plot[,snp.id],cex=1, pos=2)
        abline(h=qtl_sig_level[qtl.trait==data_qtl[,qtl.trait],sig_Z], lty=5, col='black')
        axis(2)
      }else if(var_plot == "snp_pbs.r"){
        plot(popgen_snp_plot[,.(Pos,i.pbs.r)],axes=F,pch=popgen_snp_pch, cex=1, xlim=c(data_qtl[1,qtl.focal.region.start], data_qtl[1,qtl.focal.region.end]), xlab = NULL, ylab = var_plot, col=col[i])
        text(popgen_label[,.(Pos,i.pbs.r)], labels = popgen_label[,paste(substr(GeneID,14,nchar(GeneID)),gene.name)],cex=1, pos=2)
        rect(snp.sig.gene.merge[,start], 0,  snp.sig.gene.merge[,end],5, border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.4))
        axis(2)
      }else if(var_plot == "snp_pbs.g"){
        plot(popgen_snp_plot[,.(Pos,i.pbs.g)],axes=F,pch=popgen_snp_pch, cex=1, xlim=c(data_qtl[1,qtl.focal.region.start], data_qtl[1,qtl.focal.region.end]), xlab = NULL, ylab = var_plot, col=col[i])
        text(popgen_label[,.(Pos,i.pbs.g)], labels = popgen_label[,paste(substr(GeneID,14,nchar(GeneID)),gene.name)],cex=1, pos=2)
        rect(snp.sig.gene.merge[,start], 0,  snp.sig.gene.merge[,end],5, border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.4))
        axis(2)
      }else{
        var_95<-stats95[,get(var_plot)]
        plot(dat[,.(((start+stop)/2),get(var_plot))], pch=gene_pch, cex = gene_cex, axes = F, xlim=c(data_qtl[1,qtl.focal.region.start], data_qtl[1,qtl.focal.region.end]), xlab = NULL, ylab = paste0("gene_",var_plot), col=col[i])
        abline(h=var_95, lty=5, col='black')
        axis(2)
        if(dat[avg.sig.gene==TRUE,.N]>0){
          text(dat[avg.sig.gene==TRUE,.(((start+stop)/2),get(var_plot))], labels = dat[avg.sig.gene==TRUE,paste(substr(GeneID,14,nchar(GeneID)),gene.name)],cex=1, pos=2)
          rect(avg.sig.gene.merge[,start], 0,  avg.sig.gene.merge[,end],5, border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.4))
          }
      }
      ticks<-seq(0,4e7,1e6)
      axis(1, at=ticks, labels = ticks/1e6)
    }
    mtext(sprintf("%s QTL - Chromosome %s (unit: Mb)", dat[,qtl.trait], dat[,LGn]),side=1,line=1, outer = TRUE)
  }
}  


##### Plot each QTL
{
  # create a list for quantile 95% of pbs and significance level of qtl marker.
  gene_statstoplot<-c("pbs.r","pbs.g","dp.GR")
  gene_stats95<-focal_qtl_gene[,lapply(.SD,quantile,probs=0.95,na.rm=TRUE),.SDcols = gene_statstoplot]
  qtl_sig_level<-focal_qtl_snp[sig==FALSE,.(sig_Z=max(abs_Z)),by=qtl.trait]
   
  # to plot single graph 
  plot.qtl(data_gene=focal_qtl_gene, 
           data_qtl=focal_qtl[1], 
           data_qtl_snp=focal_qtl_snp,
           data_popgen=focal_popgen_filter_plot, 
           statstoplot=c("qtl","snp_pbs.r","snp_pbs.g",gene_statstoplot), 
           stats95=gene_stats95,
           qtl_sig_level=qtl_sig_level)
  
  # plot all QTLs
  focal_qtl_plot<-focal_qtl
  for (j in 1:focal_qtl_plot[,.N]) {
    png(filename=sprintf("./Analyses_Combined/results/Foen/pbs.g_pbs.r/%s_LG%s_Start%s.png", focal_qtl_plot[j,qtl.trait],focal_qtl_plot[j,LGn], focal_qtl_plot[j,qtl.focal.region.start]),
    width = 5000, height = 3500,res=500)
    plot.qtl(data_gene=focal_qtl_gene, 
             data_qtl=focal_qtl[j], 
             data_qtl_snp=focal_qtl_snp,
             data_popgen=focal_popgen_filter_plot, 
             statstoplot=c("qtl","snp_pbs.r","snp_pbs.g",gene_statstoplot), 
             stats95=gene_stats95,
             qtl_sig_level=qtl_sig_level)
    dev.off()
  } 
}

