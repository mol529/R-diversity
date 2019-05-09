library(vegan)
##########################################################################
###########统计百分比,第一步：NTI#########################################
##########################################################################
NTI_Summ <- as.data.frame(matrix(NA,nrow = 6,ncol = 7))
colnames(NTI_Summ) <- c("Total","No_nonselection","No_homogeneous_selection","No_variable_selection","No_Dispersal_Limitation","No_Homogenizing_Dispersal","NO_Drift")

tmpn <- "02bNTI_parallel_from_pkusever/AWF_betaNTI.txt"
bNTI_AWF <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "02bNTI_parallel_from_pkusever/AWP_betaNTI.txt"
bNTI_AWP <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "02bNTI_parallel_from_pkusever/AS_betaNTI.txt"
bNTI_AS <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "02bNTI_parallel_from_pkusever/SWF_betaNTI.txt"
bNTI_SWF <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "02bNTI_parallel_from_pkusever/SWP_betaNTI.txt"
bNTI_SWP <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "02bNTI_parallel_from_pkusever/SS_betaNTI.txt"
bNTI_SS <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)


n=1
for (i in c("SWF","SWP","SS","AWF","AWP","AS")) {
  ntitmp <- as.dist(get(paste('bNTI_',i,sep = '')))
  
  No_total <- length(ntitmp)
  No_less_exp <- length(which(ntitmp < -2))
  No_more_exp <- length(which(ntitmp > 2))
  No_same_exp <- (No_total - No_less_exp - No_more_exp )

  
  NTI_Summ[n,1:4] <- c(No_total,No_same_exp,No_less_exp,No_more_exp)
  rownames(NTI_Summ)[n] <- i
  n <- n+1
}


#write.table(NTI_Summ,file = 'Summary_betaNTI.txt',quote = F,row.names = T,col.names = T,sep = "\t")




##########################################################################
###########统计百分比,第二步：RCbray######################################
##########################################################################
tmpn <- "03RCbray/SWF_RCbray.txt"
RCb_SWF <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "03RCbray/SWP_RCbray.txt"
RCb_SWP <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "99_RCbray_result_from_pkusever/SS_RCbray.txt"
RCb_SS <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)


tmpn <- "03RCbray/AWF_RCbray.txt"
RCb_AWF <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "03RCbray/AWP_RCbray.txt"
RCb_AWP <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

tmpn <- "99_RCbray_result_from_pkusever/AS_RCbray.txt"
RCb_AS <- read.table(tmpn,header = T,row.names = 1, sep='\t',check.name = F)

for (i in c("AWP","AWF","AS","SWP","SWF","SS")) {
  RCtmp <- get(paste('RCb_',i,sep = ''))#RCb_AWP
  ntitmp <- get(paste('bNTI_',i,sep = ''))#bNTI_AWP1
  
  RCtmp[ntitmp > 2] <- 0
  RCtmp[ntitmp < -2] <- 0
  
  RCtmp <- as.dist(RCtmp)
  #rc_total <- length(which(RCtmp1 != "NTI"))
  No_DL <- length(which(RCtmp > 0.95))
  No_HD <- length(which(RCtmp < -0.95))
  No_drift <- NTI_Summ[i,2] - No_DL - No_HD
  
  
  NTI_Summ[i,5:7] <- c(No_DL,No_HD,No_drift)
  
}

psumm <- NTI_Summ[,3:7]
psumm <- decostand(psumm,method = "total")
colnames(psumm) <-  c("p_homogeneous_selection","p_variable_selection","p_Dispersal_Limitation","p_Homogenizing_Dispersal","p_Drift")  

NTI_Summ <- cbind(NTI_Summ,psumm)

write.table(NTI_Summ,file = 'Quantifying_processes.txt',quote = F,row.names = T,col.names = T,sep = "\t")
