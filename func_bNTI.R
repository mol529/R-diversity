betaNTI <- function(community,cophenetic_dis,beta.reps = 999, filename="all_sample", nocore=2){
  com <- community
  cop_all <- cophenetic_dis
  
  library(picante)
  library(parallel)
  library(abind)
  
  spe_list <- colnames(com)[which(colSums(com)>0)]
  com1 <- com[,spe_list]
  cop_all1 <- cop_all[spe_list,spe_list]
  
  ##实际betaMNTD
  bmntd.weighted.ob = as.matrix(comdistnt(com1,cop_all1,abundance.weighted=T))
  write.table(bmntd.weighted.ob,file = paste(filename,"betaMNTD.txt",sep = "_"),quote = F,row.names = T,col.names = T,sep = "\t")
  ##基于的betaMNTD
  
  rand.weighted.bMNTD.comp = array(c(-999),dim=c(nrow(com1),nrow(com1),beta.reps))

  #for (rep in 1:beta.reps) {
    
  #  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(com1,taxaShuffle(cop_all1),abundance.weighted=T,exclude.conspecifics = F))

  #  print(c(date(),rep))
    
  #}
  func <- function(n,comm,cop) {
    library(picante)
    ma <- as.matrix(comdistnt(comm,taxaShuffle(cop),abundance.weighted=T,exclude.conspecifics = F))
    print(n)
    return(ma)
    }
  
  cl <- makeCluster(nocore)
  results <- parLapply(cl,1:beta.reps,func,comm <- com1,cop<-cop_all1)
  stopCluster(cl)
  rand.weighted.bMNTD.comp <- abind(results,along = 3)
  
  
  
  
  weighted.bNTI = matrix(c(NA),nrow=nrow(com1),ncol=nrow(com1))
  for (columns in 1:(nrow(com1)-1)) {
    for (rows in (columns+1):nrow(com1)) {
      
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,]
      weighted.bNTI[rows,columns] = (bmntd.weighted.ob[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      
    };
  };
 
  rownames(weighted.bNTI) = rownames(com1)
  colnames(weighted.bNTI) = rownames(com1)
  
  #write.table(weighted.bNTI,file = paste(filename,"betaNTI.txt",sep = "_"),quote = F,row.names = T,col.names = T,sep = "\t")
  return(weighted.bNTI)
}

