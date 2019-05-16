library(data.table)
setwd("/Users/pengfoen/OneDrive - University of Connecticut/Analyses_Poolseq")
#####  Load Tajima's D for robert
{
  TajD_r_filelist = list.files(path = "./Data/TajD/subsample50x_mincount2_finished", pattern="Rob.Q15", full.names = T) 
  TajD_r <- do.call(rbind,lapply(TajD_r_filelist,function(i){fread(i, sep = "\t", header = F)}))
  colnames(TajD_r) <- c("LG", "Pos", "D.nSNPs", "x.r", "D.r")
  
  TajD_r[substr(LG, 1, 2) == "MT", LGn := 0]
  TajD_r[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  TajD_r[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  TajD_r<-TajD_r[,D.r := as.numeric(D.r)]
  #TajD_r<-TajD_r[!is.na(D.r)]
  setorder(TajD_r, LGn,Pos)
  
  ###### Load Pi
  Pi_r_filelist = list.files(path = "./Data/per_chromosom_pi", pattern="Rob.Q15", full.names = T) 
  Pi_r <- do.call(rbind,lapply(Pi_r_filelist,function(i){fread(i, sep = "\t", header = F)}))
  colnames(Pi_r) <- c("LG", "Pos", "Pi.nSNPs", "x.r", "Pi.r")
  
  Pi_r[substr(LG, 1, 2) == "MT", LGn := 0]
  Pi_r[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  Pi_r[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  Pi_r<-Pi_r[,Pi.r := as.numeric(Pi.r)]
  #Pi_r<-Pi_r[!is.na(Pi.r)]
  setkey(Pi_r, LGn, Pos)

}

{
  #####  Load Tajima's D for gos
  TajD_g_filelist = list.files(path = "./Data/TajD/subsample50x_mincount2_finished", pattern="Gos.Q15", full.names = T) 
  TajD_g <- do.call(rbind,lapply(TajD_g_filelist,function(i){fread(i, sep = "\t", header = F)}))
  colnames(TajD_g) <- c("LG", "Pos", "D.nSNPs", "x.g", "D.g")
  
  TajD_g[substr(LG, 1, 2) == "MT", LGn := 0]
  TajD_g[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  TajD_g[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  TajD_g<-TajD_g[,D.g := as.numeric(D.g)]
  #TajD_g<-TajD_g[!is.na(D.r)]
  setorder(TajD_g, LGn,Pos)
  
  
  ###### Load Pi
  Pi_g_filelist = list.files(path = "./Data/per_chromosom_pi", pattern="Gos.Q15", full.names = T) 
  Pi_g <- do.call(rbind,lapply(Pi_g_filelist,function(i){fread(i, sep = "\t", header = F)}))
  colnames(Pi_g) <- c("LG", "Pos", "Pi.nSNPs", "x.g", "Pi.g")
  
  Pi_g[substr(LG, 1, 2) == "MT", LGn := 0]
  Pi_g[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  Pi_g[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  Pi_g<-Pi_g[,Pi.g := as.numeric(Pi.g)]
  #Pi_g<-Pi_g[!is.na(Pi.r)]
  setkey(Pi_g, LGn, Pos)
}

{
  #####  Load Tajima's D for say
  TajD_s_filelist = list.files(path = "./Data/TajD/subsample50x_mincount2_finished", pattern="Say.Q15", full.names = T) 
  TajD_s <- do.call(rbind,lapply(TajD_s_filelist,function(i){fread(i, sep = "\t", header = F)}))
  colnames(TajD_s) <- c("LG", "Pos", "D.nSNPs", "x.s", "D.s")
  
  TajD_s[substr(LG, 1, 2) == "MT", LGn := 0]
  TajD_s[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  TajD_s[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  TajD_s<-TajD_s[,D.s := as.numeric(D.s)]
  #TajD_s<-TajD_s[!is.na(D.r)]
  setorder(TajD_s, LGn,Pos)
  
  
  ###### Load Pi
  Pi_s_filelist = list.files(path = "./Data/per_chromosom_pi", pattern="Gos.Q15", full.names = T) 
  Pi_s <- do.call(rbind,lapply(Pi_s_filelist,function(i){fread(i, sep = "\t", header = F)}))
  colnames(Pi_s) <- c("LG", "Pos", "Pi.nSNPs", "x.s", "Pi.s")
  
  Pi_s[substr(LG, 1, 2) == "MT", LGn := 0]
  Pi_s[substr(LG, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(LG, 6, nchar(LG))))]
  Pi_s[substr(LG, 1, 2) == "sc", LGn := as.numeric(substr(LG, 10, nchar(LG)))]
  Pi_s<-Pi_s[,Pi.s := as.numeric(Pi.s)]
  #Pi_s<-Pi_s[!is.na(Pi.r)]
  setkey(Pi_s, LGn, Pos)
}


dat<-Pi_g
var_plot<-"Pi.g"
Chr<-12

maxpos<-max(dat[LGn==Chr,Pos])
plot(dat[LGn==Chr,.(Pos/1e6,get(var_plot))],pch = 16, cex = 0.6, axes=F, xlab = paste0("Chr ",Chr), ylab = var_plot)
ticks<-seq(0,maxpos/1e6,1)
axis(1, at=seq(0,maxpos/1e6,0.5), labels = F)
axis(1, at=seq(0,maxpos/1e6,1), labels = seq(0,maxpos/1e6,1))
axis(2)