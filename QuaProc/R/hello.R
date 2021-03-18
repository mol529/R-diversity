# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

###########################################################################################################
###cophenetic for big tree
###########################################################################################################

##############################################################
####以下代码需要在服务器跑，个人电脑可能带不动！！！！########
##############################################################

bigtree_cophenetic <- function(rooted_tree, subsize = 9000){
  t0=Sys.time()
  library(picante)
  #tree <- read.tree(file = "rooted-tree.nwk")
  tree <- rooted_tree
  
  
  tips_list <- tree$tip.label
  ##设置每个子集大小,要小于2万，即sub_size小于10000。
  sub_size <- subsize
  n <- ceiling(length(tips_list)/sub_size)
  if (n < 3) {
    cop_all <- cophenetic(tree)
  }
  if (n > 2) {
    sub_tiplist <- list()
    for (i in 1:(n-1)) {
      sub_tmp <- tips_list[(1+(i-1)*sub_size) : (i*sub_size)]
      sub_tiplist[[i]] <- sub_tmp
    }
    sub_tiplist[[n]] <- tips_list[(1+(n-1)*sub_size) : length(tips_list)]
    
    cop_all <- matrix(0,nrow = length(tips_list),ncol = length(tips_list))
    cop_all <- as.data.frame(cop_all)
    rownames(cop_all) <- tips_list
    colnames(cop_all) <- tips_list
    
    for (i in 1:n) {
      for (j in i:n) {
        if (j>i) {
          a <- sub_tiplist[[i]]
          b <- sub_tiplist[[j]]
          ab <- c(a,b)
          print(paste(i,j))
          print(length(ab))
          c <- setdiff(tips_list ,ab) 
          tree1<-drop.tip(tree,c)
          cop <- cophenetic(tree1)
          cop_all[rownames(cop),colnames(cop)] <- cop[rownames(cop),colnames(cop)]
        }
        
      }
      
    }
    
  }
  
  print(paste0("Big tree distance finished time: ",Sys.time()))
  print(Sys.time()-t0)
  return(cop_all)
}




###########################################################################################################
###betaNTI
###########################################################################################################
betaNTI <- function(com,
                    cophenetic_dis,
                    beta.reps = 999, 
                    filename="all_sample", 
                    nocore=2, 
                    nbin=5,
                    mth_snowfall=F){
  print(paste0("Process begin: ",Sys.time()))
  print(paste0("This com table is a ",nrow(com)," x ",ncol(com)," matrix."))
  t_start=Sys.time()
  
  if(mth_snowfall){
    library(snowfall)
    print("Use Snowfall package")
  } else {
    library(parallel)
    print("Use Parallel package")
  }
  
  library(picante)
  #library(abind)
  library(bigmemory)
  #library(tidyr)
  
  com <- com[,colSums(com)>0]
  cophenetic_dis <- cophenetic_dis[colnames(com),colnames(com)]
  
  ##实际betaMNTD
  print(paste0("Begin betaMNTD: ",Sys.time()))
  t0=Sys.time()
  bmntd.weighted.ob = picante::comdistnt(com,cophenetic_dis,abundance.weighted=T)
  bmntd.weighted.ob = as.matrix(bmntd.weighted.ob)
  bmntd.weighted.ob[upper.tri(bmntd.weighted.ob)] <- NA
  diag(bmntd.weighted.ob) <- NA
  bmntd.weighted.ob <- as.data.frame(bmntd.weighted.ob)
  bmntd.weighted.ob$id2 <- rownames(bmntd.weighted.ob)
  bmntd.weighted.ob <- tidyr::gather(bmntd.weighted.ob,key = "id1", value = "disval",na.rm = T,-id2)
  write.table(bmntd.weighted.ob,file = paste(filename,"betaMNTD.txt",sep = "_"),quote = F,row.names = F,col.names = T,sep = "\t")
  print(paste0("betaMNTD time:",format(Sys.time()-t0)))
  
  ##基于的betaMNTD
  
  #rand.weighted.bMNTD.comp = array(c(-999),dim=c(nrow(com),nrow(com),beta.reps))
  
  #for (rep in 1:beta.reps) {
  
  #  rand.weighted.bMNTD.comp[,,rep] = as.matrix(picante::comdistnt(com,taxaShuffle(cophenetic_dis),abundance.weighted=T,exclude.conspecifics = F))
  
  #  print(c(date(),rep))
  
  #}
  func <- function(n,comm,cop) {
    library(picante)
    #library(tidyr)
    ma <- as.matrix(picante::comdistnt(comm,taxaShuffle(cop),abundance.weighted=T,exclude.conspecifics = F))
    ma[upper.tri(ma)] <- NA
    diag(ma) <- NA
    #print(n)
    colnames(ma) <- rownames(ma) <- 1:nrow(ma)
    ma <- as.data.frame(ma)
    ma$id2 <- rownames(ma)
    ma <- tidyr::gather(ma,key = "id1", value = "disval",na.rm = T,-id2)
    gc()
    return(ma)
  }
  
  print(paste0("Begin null model: ",Sys.time()))
  nr=nrow(com)*(nrow(com)-1)/2
  
  if(mth_snowfall){
    sfInit(parallel = TRUE, cpus = nocore)
    sfLibrary(picante)
    sfLibrary(tidyr)
    sfExport("com","cophenetic_dis")
  } else {
    cl <- makeCluster(nocore)
  }
  
  if(nbin == 1){
    t0=Sys.time()
    if(mth_snowfall){
      results <- sfLapply(1:beta.reps,func,comm <- com,cop<-cophenetic_dis) 
    } else {
      results <- parLapply(cl,1:beta.reps,func,comm <- com,cop<-cophenetic_dis)
    }
    
    rand.weighted.bMNTD.comp <- abind::abind(results,along = 1)
    rm(results)
    gc()
    print(paste0("null model time:",format(Sys.time()-t0)))
  } else {
    subn <- round(beta.reps/nbin,digits = 0)
    rand.weighted.bMNTD.comp <- big.matrix(nrow=beta.reps*nr, ncol=3, init=-999, backingfile='parall_null_model_res.bin', descriptorfile='parall_null_model_res.desc')
    
    for(i in 1:nbin){
      print(paste0("start bin",i))
      t0=Sys.time()
      if(i != nbin){
        print(paste0((i-1)*subn+1,"-",(i*subn)))
        if(mth_snowfall){
          results <- sfLapply(((i-1)*subn+1):(i*subn),func,comm <- com,cop<-cophenetic_dis)
        } else {
          results <- parLapply(cl,((i-1)*subn+1):(i*subn),func,comm <- com,cop<-cophenetic_dis)
        }
        rand.weighted.bMNTD.comp[(nr*(i-1)*subn+1):(nr*i*subn),] <- abind::abind(results,along = 1) 
      } else {
        print(paste0((i-1)*subn+1,"-",beta.reps))
        if(mth_snowfall){
          results <- sfLapply(((i-1)*subn+1):beta.reps,func,comm <- com,cop<-cophenetic_dis)
        } else {
          results <- parLapply(cl,((i-1)*subn+1):beta.reps,func,comm <- com,cop<-cophenetic_dis)
        }
        rand.weighted.bMNTD.comp[(nr*(i-1)*subn+1):(nr*beta.reps),] <- abind::abind(results,along = 1)
      } 
      print(paste0("bin",i,"null model time:",format(Sys.time()-t0)))
      rm(results)
      gc()
    }
  } 
  
  if(mth_snowfall){
    sfStop()
  } else {
    stopCluster(cl)
  }
  
  print("Parallel process fishned")
  
  rm(cophenetic_dis)
  gc()
  
  print(paste0("Begin loops: ",Sys.time()))
  weighted.bNTI = bmntd.weighted.ob
  weighted.bNTI[,3] <- -999
  t0=Sys.time()
  for (id1 in 1:(nrow(com)-1)) {
    for (id2 in (id1+1):nrow(com)) {
      rand.vals = as.numeric(rand.weighted.bMNTD.comp[rand.weighted.bMNTD.comp[,2]==id1 & rand.weighted.bMNTD.comp[,1]==id2,3])
      
      weighted.bNTI[weighted.bNTI[,"id1"]==rownames(com)[id1] & weighted.bNTI[,"id2"]==rownames(com)[id2],3] <- 
        (bmntd.weighted.ob[weighted.bNTI[,"id1"]==rownames(com)[id1] & weighted.bNTI[,"id2"]==rownames(com)[id2],3] - 
           mean(rand.vals)) / sd(rand.vals)
    }
  }
  
  colnames(weighted.bNTI) <- c("Sample1","Sample2","betaNTI")
  print(paste0("loops  time:",format(Sys.time()-t0)))
  write.table(weighted.bNTI,file = paste(filename,"betaNTI.txt",sep = "_"),quote = F,row.names = F,col.names = T,sep = "\t")
  #return(weighted.bNTI)
  print(paste0("Total time:",format(Sys.time()-t_start)))
  
}


###########################################################################################################
###RCbray
###########################################################################################################
raup_crick_modified=function(comm,
                             nocore=2,
                             reps=9999,
                             filename="all_sample",
                             nbin=5,
                             mth_snowfall=F){
  
  ##comm: a species by site matrix, with row names for plots.  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model
  print(paste0("RCbray begin: ",Sys.time()))
  print(paste0("This com table is a ",nrow(comm)," x ",ncol(comm)," matrix."))
  t_start = Sys.time()
  
  if(mth_snowfall){
    library(snowfall)
    print("Use Snowfall package")
  } else {
    library(parallel)
    print("Use Parallel package")
  }
  
  #library(vegan)
  #library(tidyr)
  library(bigmemory)
  
  #对物种求和
  com_sum <- colSums(comm)
  
  ## count number of sites and total species richness across all plots，提取物种名称 (gamma)
  n_sites <-nrow(comm)
  n_sp <- ncol(comm)
  gamma <- colnames(comm)
  
  ##make the comm matrix into a pres/abs. (overwrites initial comm matrix):
  comm_standard <- vegan::decostand(comm,"pa")
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- colSums(comm_standard)
  
  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(rowSums(comm_standard)))
  
  #generate the matrix for the input of func_subRC function
  mtr <- matrix(0,length(alpha_levels),length(alpha_levels))
  rownames(mtr) <- colnames(mtr) <- 1:length(alpha_levels)
  mtr[upper.tri(mtr)] <- NA
  mtr <- as.data.frame(mtr)
  mtr$id2 <- rownames(mtr)
  long_bray <- tidyr::gather(mtr,key = "id1", value = "repsno",na.rm = T,-id2)
  
  #a1 <- rep(long_bray$id1,reps)
  #a2 <- rep(long_bray$id2,reps)
  #repsno <- rep(1:reps, each = nrow(long_bray))
  #vfsmmy <- data.frame(a1,a2,repsno)
  vfsmmy <- matrix(c(rep(long_bray$id1,reps),rep(long_bray$id2,reps),rep(1:reps, each = nrow(long_bray))),reps*nrow(long_bray),3)
  rm(long_bray,mtr)
  gc()
  
  #构建用于并行函数
  #vf: containing 3 factors needed for cycle
  func_subRC <- function(nl,vfsmmy,comm,alpha_levels,occur){
    #library(snowfall)
    #library(vegan)	 
    a1 <- as.numeric(vfsmmy[nl,1])
    a2 <- as.numeric(vfsmmy[nl,2])
    #n <- vfsmmy[nl,3]
    
    n_sp <- ncol(comm)
    gamma <- colnames(comm)
    com_sum <- colSums(comm)
    ##two empty null communities of size gamma:
    com1 <- com2 <- rep(0,n_sp)
    names(com1) <- names(com2) <- gamma
    
    ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
    com1[sample(gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
    ##same for com2:
    com2[sample(gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
    
    ##how many species are shared in common?
    #null_dist[i]<-sum((com1+com2)>1)
    
    ##有两种方法，测试证明，两者结果相似，mantel test结果R为99.1%，p为0.001.这里选择方法2。
    ###方法一：用所有样品平均多度作为概率抽取物种多度（会出现一部分物种多度为0）
    
    #allreads1 <- sum(comm[a1,])
    #com1list <- which(com1==1)
    #com1_pro <- com_sum[com1list]
    #ind1 <- as.factor(sample(com1list,allreads1,replace=TRUE, prob=com1_pro))
    #ind1_summary <- t(as.matrix(table(ind1)))
    
    #allreads2 <- sum(comm[a2,])
    #com2list <- which(com2==1)
    #com2_pro <- com_sum[com2list]
    #ind2 <- as.factor(sample(com2list,allreads2,replace=TRUE, prob=com2_pro))
    #ind2_summary <- t(as.matrix(table(ind2)))
    
    
    #pairma <- matrix(0,nrow=2,ncol=n_sp)
    #colnames(pairma) <- 1:n_sp
    #pairma[1,colnames(ind1_summary)] <- ind1_summary
    #pairma[2,colnames(ind2_summary)] <- ind2_summary
    
    #dis1 <- vegan::vegdist(pairma,method="bray")
    
    ###方法二：以所有样品平均多度为比例，计算物种多度（保证了所有样品多度不为0）
    allreadcount <- sum(comm[1,])
    
    com1list <- which(com1==1)
    com1_pro <- com_sum[com1list]
    com1_pro <- allreadcount*com1_pro/sum(com1_pro)
    
    com2list <- which(com2==1)
    com2_pro <- com_sum[com2list]
    com2_pro <- allreadcount*com2_pro/sum(com2_pro)
    
    pairma <- matrix(0,nrow=2,ncol=n_sp)
    colnames(pairma) <- gamma
    pairma[1,names(com1_pro)] <- com1_pro
    pairma[2,names(com2_pro)] <- com2_pro
    
    dis2 <- vegan::vegdist(pairma,method="bray")
    
    output <- c(alpha_levels[a1],alpha_levels[a2],dis2)
    
    return(output)
  }
  
  
  ##make_null:
  print(paste0("RCbray null model begin: ",Sys.time()))
  
  if(mth_snowfall){
    sfInit(parallel = TRUE, cpus = nocore)
    sfLibrary(vegan)
    sfExport("vfsmmy","comm", "alpha_levels", "occur")
  } else {
    cl <- makeCluster(nocore)
  }
  
  ttl_reps <- nrow(vfsmmy)  
  
  
  if(nbin == 1) {
    t0=Sys.time()
    
    if(mth_snowfall){
      results <- sfLapply(1:ttl_reps,func_subRC,vfsmmy <- vfsmmy, comm <- comm,alpha_levels<-alpha_levels,occur <- occur)
    } else {
      results <- parLapply(cl,1:ttl_reps,func_subRC,vfsmmy <- vfsmmy, comm <- comm,alpha_levels<-alpha_levels,occur <- occur)
    }
    
    null_dist <- matrix(unlist(results), nrow=3)

    rm(results)
    gc()
    print(paste0("null model time:",format(Sys.time()-t0,digits=4)))
  } else {
    subn <- round(ttl_reps/nbin,digits = 0)
    print(paste0("Need ",nbin, " loops. Each loop has ",subn," reps." ))
    null_dist <- big.matrix(nrow=3, ncol=ttl_reps, init=-999, backingfile='parall_null_model_res.bin', descriptorfile='parall_null_model_res.desc')

    
    for(i in 1:nbin){
      print(paste0("start bin",i))
      t0=Sys.time()
      if(i != nbin){
        if(mth_snowfall){
          results <- sfLapply(((i-1)*subn+1):(i*subn),func_subRC,vfsmmy <- vfsmmy, comm <- comm,alpha_levels<-alpha_levels,occur <- occur)
        } else {
          results <- parLapply(cl,((i-1)*subn+1):(i*subn),func_subRC,vfsmmy <- vfsmmy, comm <- comm,alpha_levels<-alpha_levels,occur <- occur)
        }
        
        null_dist[,((i-1)*subn+1):(i*subn)] <- matrix(unlist(results), nrow=3)
      } else {
        if(mth_snowfall){
          results <- sfLapply(((i-1)*subn+1):ttl_reps,func_subRC,vfsmmy <- vfsmmy, comm <- comm,alpha_levels<-alpha_levels,occur <- occur)
        } else {
          results <- parLapply(cl,((i-1)*subn+1):ttl_reps,func_subRC,vfsmmy <- vfsmmy, comm <- comm,alpha_levels<-alpha_levels,occur <- occur)
        }
        
        null_dist[,((i-1)*subn+1):ttl_reps] <- matrix(unlist(results), nrow=3)
      } 
      print(paste0("bin",i,"null model time:",format(Sys.time()-t0,digits=4)))
      rm(results)
      gc()
    }
  }
  
  #sfStop()
  
  print("Parallel null model processes finished")
  
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  
  
  disbray <- as.matrix(vegan::vegdist(comm,method="bray"))
  disbray[upper.tri(disbray)] <- NA
  diag(disbray) <- NA
  disbray <- as.data.frame(disbray)
  disbray$id1 <- rownames(disbray)
  disbray <- tidyr::gather(disbray,key = "id2", value = "rn",na.rm = T,-id1)
  
  rm(occur, comm,vfsmmy)
  gc()
  
  func_mtr <- function(nlm, disbray, reps, comm_standard,null_dist){
    #require(bigmemory)
    ##how many species are shared between the two sites:
    n_shared_obs <- disbray[nlm,3]
    
    ## what was the observed richness of each site?
    obs_a1<-sum(comm_standard[disbray[nlm,1],])
    obs_a2<-sum(comm_standard[disbray[nlm,2],])
    
    ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
    obs_a_pair<-sort(c(obs_a1, obs_a2))
    
    nullrc <- null_dist[3,null_dist[1,]==obs_a_pair[1] & null_dist[2,]==obs_a_pair[2]]
    
    ##how many null observations is the observed value tied with?
    num_exact_matching_in_null<-sum(nullrc==n_shared_obs)
    
    ##how many null values are bigger than the observed value?
    num_greater_in_null<-sum(nullrc>n_shared_obs)
    
    rc <- ((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)			
    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
    rc <- (rc - 0.5)*2
    ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
    rc<- -1*rc
    
    return(c(disbray[nlm,1],disbray[nlm,2],round(rc, digits=3)))
  }
  
  
  null_dist <- as.matrix(null_dist)
  
  if(mth_snowfall){
    #sfInit(parallel = TRUE, cpus = nocore)
    sfExport("disbray","reps", "comm_standard", "null_dist")
    rcres <- sfLapply(1:nrow(disbray),func_mtr,
                      disbray <- disbray,
                      reps <- reps,
                      comm_standard<-comm_standard,
                      null_dist <- null_dist)
    sfStop()
  } else {
    rcres <- parLapply(cl,1:nrow(disbray),func_mtr,
                       disbray <- disbray,
                       reps <- reps,
                       comm_standard<-comm_standard,
                       null_dist <- null_dist)
    stopCluster(cl)
  }
  
  longRC <- t(matrix(unlist(rcres), nrow=3))
  colnames(longRC) <- c("Sample1","Sample2","RCbray")
  
  write.table(longRC,file = paste(filename,"RCbray.txt",sep = "_"),quote = F,row.names = F,col.names = T,sep = "\t")
  print(paste0("Total time:",format(Sys.time()-t_start,,digits=4)))
}





