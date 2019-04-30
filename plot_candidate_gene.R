#################################################################
#  1.  Load packages
#################################################################
{
  library(data.table)
}

#################################################################
#  2.  Load data
#################################################################
rm(list=ls())

## function to convert linkage group names from string to numeric number
convert_LG_name<-function(dataset, LG_column){
  dataset[substr(get(LG_column), 1, 2) == "MT", LGn := 0]
  dataset[substr(get(LG_column), 1, 2) == "gr", LGn := as.numeric(as.roman(substr(get(LG_column), 6, nchar(get(LG_column)))))]
  dataset[substr(get(LG_column), 1, 2) == "sc", LGn := as.numeric(substr(get(LG_column), 10, nchar(get(LG_column))))]
}

{
  #setwd("D:/Foen Peng/OneDrive - University of Connecticut/")
  setwd("/Users/pengfoen/OneDrive - University of Connecticut/")
  
  gene_loc<-unique(fread("./Analysis_Expression/Data files RAW/GeneLocations.csv",sep=","))
  gene_name<-fread("./Analysis_Expression/Data files RAW/GeneID_to_GeneName.csv")
  gene_LG<-fread("./Analysis_Expression/Data files RAW/GeneID to LinkageGroup.csv",header = TRUE)
  setnames(gene_name,c("Ensembl Gene ID","Associated Gene Name"),c("GeneID","gene.name"))
  convert_LG_name(gene_LG, "LinkageGroup")
  # LG file has a space at the begining of GeneID.
  gene_LG[,GeneID:=substring(V1,2)]
  
  setkey(gene_loc,GeneID)
  setkey(gene_LG, GeneID)
  setkey(gene_name,GeneID)
  
  # merge three tables to gene_info
  gene_info<-gene_LG[gene_loc, on=c("GeneID")]
  gene_info<-gene_info[gene_name]
  gene_info<-gene_info[,!c("V1", "LinkageGroup")]
  setkey(gene_info,LGn,Start,Stop)
  
  # seperate the full length from the exons
  gene_info[,full:=FALSE]
  gene_info[,c("gene_start","gene_end"):=list(min(Start),max(Stop)),by=GeneID]
  gene_info[Start==gene_start&Stop==gene_end,full:=TRUE]
  gene_info[,c("gene_start","gene_end"):=NULL]
  
  region_to_plot<-fread("./Analyses_Combined/results/Foen/region_roomin.csv")
  gene_to_plot<-fread("./Analyses_Combined/results/Foen/candidate_list_all_evidence.csv")[priority_level==1,GeneID]
  popgen_all<-fread("./Analyses_Poolseq/popgen_info_filter_seq_depth.csv",header = T)
  
}

#################################################################
#  3.  plot gene
#################################################################

plot_gene <- function(dat, focalLG, specifybps = F, minpos, maxpos, gene.name, gene.id=gene.id, exon,snp.exon, statstoplot){
  npops <- length(statstoplot)
  par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))

  options(scipen=5)
  col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
  for(i in 1:length(statstoplot)){
    if(statstoplot[i]=="dp.GR"){
      y_lim=c(0,1.1)
    }else{
      y_lim=c(0,5.5)
    }
    
    plot(dat[,.(Pos,eval(parse(text = statstoplot[i])))], pch = 16, cex = 1,type = 'o', axes = T, ylim=y_lim, ylab = statstoplot[i],col=col[i])
    rect(exon[,Start], 0,  exon[,Stop],y_lim[2]*(1/1.1), border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.2))
    if(length(snp.exon)!=0){
      points(snp.exon,rep(y_lim[2]*(1/1.1),length(snp.exon)))
      text(snp.exon,y_lim[2]*(1/1.1),labels=substr(snp.exon, nchar(snp.exon)-4+1, nchar(snp.exon)),cex=0.7, pos=3)
    }

    #abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),((stop-start)/2)/1000000],col=rgb(red=0,green=0,blue=0,alpha=100,maxColorValue = 255))
  }
  mtext(paste("Candidate",gene.name,gene.id),side=1,line=1, outer = TRUE)
}

gene_to_plot<-c("ENSGACG00000015525","ENSGACG00000015530","ENSGACG00000015539")

for (i in gene_to_plot) {
  gene_info_plot<-gene_info[GeneID==i,]
  print(i)
  gene_start_plot<-gene_info_plot[full==TRUE,Start] -3000
  gene_end_plot<-gene_info_plot[full==TRUE,Stop] + 3000
  gene_name_plot<-gene_info_plot[full==TRUE,gene.name]
  LGn_plot<-gene_info_plot[full==TRUE,LGn]
  exon_plot<-gene_info_plot[,.(Start,Stop)]
  exon_plot[,length:=Stop-Start]
  if(exon_plot[,.N]==1){
    exon_plot[,EXON:=T]
  }else{
    exon_plot[length<0.8*max(length), EXON:=T]
  }
  SNP_plot<-popgen_all[LGn==LGn_plot & Pos >= gene_start_plot & Pos <= gene_end_plot,.SD,keyby=Pos]
  SNP_exon<-SNP_plot[Pos %inrange% exon_plot[EXON==T, .(Start,Stop)]][dp.GR>0.8]
  print(SNP_exon[,.(Pos, dp.GR, i.Base_A)])
  #SNP_Exon<-exon_plot[EXON==T][SNP_plot,on=c("Start<Pos","Stop>Pos"),nomatch=0L, allow.cartesian=F][dp.GR>0.8,]
  
  #png(filename=sprintf("./Analyses_Combined/results/Foen/Candidate_gene/%s_%s.png",substr(i,14,nchar(i)),gene_name_plot),width = 3000, height = 2000,res=300)
  #plot_gene(SNP_plot,focalLG=LGn_plot, minpos=gene_start_plot, maxpos=gene_end_plot, gene.name=gene_name_plot,gene.id=i, exon=exon_plot[EXON==T],snp.exon=SNP_exon[,Pos],
  #         statstoplot = c("dp.GR","pbs.r","pbs.g"))
  #dev.off()
} 

#################################################################
#  4.  plot region of interests
#################################################################
plot_region <- function(dat, focalLG, minpos, maxpos, gene.name, gene.id, exon, statstoplot){
  npops <- length(statstoplot)
  par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
  
  options(scipen=5)
  col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
  for(i in 1:length(statstoplot)){
    if(statstoplot[i] %in% c("dp.GR","GRFst.nonneg")){
      y_lim=c(0,1)
    }else{
      y_lim=c(0,5)
    }
    plot(dat[,.(Pos,eval(parse(text = statstoplot[i])))], pch = 16, cex = 0.5, axes = F, ylim=y_lim, ylab = statstoplot[i],col=col[i])
    rect(exon[,Start], 0,  exon[,Stop],y_lim[2], border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.1))
    text(exon[,(Start+Stop)/2],y_lim[2],labels = paste(substr(gene.id,14,nchar(gene.id)),gene.name,sep="\n"),cex=0.7, pos=1)
    axis(1, at=seq(minpos,maxpos,0.025E6), labels = F)
    axis(1, at=seq(minpos,maxpos,0.05E6), labels = seq(minpos/1e6,maxpos/1e6,0.05))
    axis(2)
  }
  mtext(paste("LG",focalLG,":", minpos,"-",maxpos),side=1,line=1, outer = TRUE)
}

region_to_plot<-data.table(
  LGn = 2,
  start = 7500000,
  end = 8000000
)
for (i in 1:region_to_plot[,.N]) {
  LGn_plot<-region_to_plot[i,LGn]
  region_start_plot<-region_to_plot[i,start]
  region_end_plot<-region_to_plot[i,end]
  
  png(filename=sprintf("./Analyses_Poolseq/Results/Foen/Region_zoomin/LG%s_%s-%s.png",LGn_plot,region_start_plot/1e6,region_end_plot/1e6),width = 3000, height = 2000,res=300)

  gene_info_plot<-gene_info[Start>region_start_plot & Stop<region_end_plot & LGn==LGn_plot,]
  gene_name_plot<-gene_info_plot[full==TRUE,gene.name]
  gene_id_plot<-gene_info_plot[full==TRUE,GeneID]
  exon_plot<-gene_info_plot[full==TRUE,.(Start, Stop)]
  SNP_plot<-popgen_all[LGn==LGn_plot & Pos >= region_start_plot & Pos <= region_end_plot,.SD,keyby=Pos]
  plot_region(SNP_plot,focalLG=LGn_plot, minpos=region_start_plot, maxpos=region_end_plot, gene.name=gene_name_plot,gene.id=gene_id_plot, exon=exon_plot,
            statstoplot = c("dp.GR","GRFst.nonneg", "pbs.r","pbs.g"))
  dev.off()
} 
