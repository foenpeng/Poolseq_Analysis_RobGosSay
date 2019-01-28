#################################################################
#  1.  Load packages
#################################################################
{
  library(data.table)
}

#################################################################
#  2.  Load data
#################################################################
{
  setwd("/Users/pengfoen/OneDrive - University of Connecticut/Analyses_QTL")
  qtl_fib<-fread("./result/SingleMarkerResultsDan/Fibrosis_SingleMarkerResults.csv", header=TRUE)
  qtl_snp_pos<- fread("./result/SingleMarkerResultsDan/snp_pos_QTL.csv", header=TRUE)
  
  
  qtl_snp_pos[substr(V1, 1, 2) == "gr", LGn := as.numeric(as.roman(substr(V1, 6, nchar(V1))))]
}
