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
bigtree_cophenetic <- function(rooted_tree, subsize = 9000){
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

  return(cop_all)
}



###########################################################################################################
###betaNTI
###########################################################################################################
betaNTI <- function(community,cophenetic_dis,beta.reps = 999, filename="all_sample", nocore=2,returnresult=F){
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

  write.table(weighted.bNTI,file = paste(filename,"betaNTI.txt",sep = "_"),quote = F,row.names = T,col.names = T,sep = "\t")

  if (returnresult == T){
  return(weighted.bNTI)
    }
}




###########################################################################################################
###RCbray
###########################################################################################################
raup_crick_modified=function(spXsite, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=9999, nocore=2, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){

  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.


  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model


  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  library(vegan)
  library(parallel)

  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }

  #对物种求和
  com_sum <- apply(spXsite,2,sum)

  ## count number of sites and total species richness across all plots，提取物种名称 (gamma)
  n_sites<-nrow(spXsite)
  n_sp <- ncol(spXsite)
  gamma <- colnames(spXsite)


  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite_standard

  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite_standard, MARGIN=2, FUN=sum)


  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur<-rep(1,n_sp)
  }





  #构建并行函数
  func_subRC <- function(n,spXsite,alpha_levels,occur,a1,a2){
    library(parallel)
    library(vegan)
    n_sp <- ncol(spXsite)
    gamma <- colnames(spXsite)
    com_sum <- apply(spXsite,2,sum)
    ##two empty null communities of size gamma:
    com1<-rep(0,n_sp)
    com2<-rep(0,n_sp)
    names(com1) <- gamma
    names(com2) <- gamma
    ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
    com1[sample(gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1

    ##same for com2:
    com2[sample(gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1

    ##how many species are shared in common?
    #null_shared_spp[i]<-sum((com1+com2)>1)

    ##有两种方法，测试证明，两者结果相似，mantel test结果R为99.1%，p为0.001.这里选择方法2。
    ###方法一：用所有样品平均多度作为概率抽取物种多度（会出现一部分物种多度为0）

    #allreads1 <- sum(spXsite[a1,])
    #com1list <- which(com1==1)
    #com1_pro <- com_sum[com1list]
    #ind1 <- as.factor(sample(com1list,allreads1,replace=TRUE, prob=com1_pro))
    #ind1_summary <- t(as.matrix(table(ind1)))

    #allreads2 <- sum(spXsite[a2,])
    #com2list <- which(com2==1)
    #com2_pro <- com_sum[com2list]
    #ind2 <- as.factor(sample(com2list,allreads2,replace=TRUE, prob=com2_pro))
    #ind2_summary <- t(as.matrix(table(ind2)))


    #pairma <- matrix(0,nrow=2,ncol=n_sp)
    #colnames(pairma) <- 1:n_sp
    #pairma[1,colnames(ind1_summary)] <- ind1_summary
    #pairma[2,colnames(ind2_summary)] <- ind2_summary

    #dis1 <- vegdist(pairma,method="bray")

    ###方法二：以所有样品平均多度为比例，计算物种多度（保证了所有样品多度不为0）
    allreads1 <- sum(spXsite[a1,])
    com1list <- which(com1==1)
    com1_pro <- com_sum[com1list]
    com1_pro <- allreads1*(com1_pro/sum(com1_pro))

    allreads2 <- sum(spXsite[a2,])
    com2list <- which(com2==1)
    com2_pro <- com_sum[com2list]
    com2_pro <- allreads2*(com2_pro/sum(com2_pro))


    pairma <- matrix(0,nrow=2,ncol=n_sp)
    colnames(pairma) <- gamma
    pairma[1,names(com1_pro)] <- com1_pro
    pairma[2,names(com2_pro)] <- com2_pro

    dis2 <- vegdist(pairma,method="bray")


    return(dis2)
  }



  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(apply(spXsite_standard, MARGIN=1, FUN=sum)))

  ##make_null:

  ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels.

  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1

  ##null_array will hold the actual null distribution values.  Each element of the array corresponds to a null distribution for each combination of alpha values.  The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  Later the function will find the row of alpha_table with the right combination of alpha values.  That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels.
  null_array<-list()

  ##looping over each combination of alpha levels:
  for(b1 in 1:length(alpha_levels)){
    for(b2 in b1:length(alpha_levels)){

      cl <- makeCluster(nocore)
      results <- parLapply(cl,1:reps,func_subRC,spXsite <- spXsite,alpha_levels<-alpha_levels,occur <- occur,a1 <- b1,a2 <- b2)
      stopCluster(cl)
      null_shared_spp <- as.matrix(rbind(results))


      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp

      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[b1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[b2]

      #increment the counter for the columns of the alpha table/ elements of the null array
      print(paste(b1,b2,sep = '_'))

      col_count<-col_count+1

    }
  }

  ##create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")


  #####################
  ##do the test:



  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))


  ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){

      ##how many species are shared between the two sites:
      #n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      n_shared_obs <- vegdist(spXsite[c(i,j),],method="bray")

      ## what was the observed richness of each site?
      obs_a1<-sum(spXsite_standard[i,])
      obs_a2<-sum(spXsite_standard[j,])

      ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))

      ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))

      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)

      ##how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)



      rc<-(num_greater_in_null)/reps




      if(split_ties){

        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }



      if(!classic_metric){

        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1

        rc<-(rc-.5)*2
      }


      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }

      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }


      ##store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)


    }
  }


  if(as.distance.matrix){
    results<-as.dist(results)
  }


  return(results)

}


#######################--------------------------
#RC_split
#--------------------------------------------
raup_crick_modified_split=function(spXsite,startn,endn, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=9999, nocore=2, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){

  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.


  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model


  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  library(vegan)
  library(parallel)

  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }

  #对物种求和
  com_sum <- apply(spXsite,2,sum)

  ## count number of sites and total species richness across all plots，提取物种名称 (gamma)
  n_sites<-nrow(spXsite)
  n_sp <- ncol(spXsite)
  gamma <- colnames(spXsite)


  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite_standard

  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite_standard, MARGIN=2, FUN=sum)


  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur<-rep(1,n_sp)
  }





  #构建并行函数
  func_subRC <- function(n,spXsite,alpha_levels,occur,a1,a2){
    library(parallel)
    library(vegan)
    n_sp <- ncol(spXsite)
    gamma <- colnames(spXsite)
    com_sum <- apply(spXsite,2,sum)
    ##two empty null communities of size gamma:
    com1<-rep(0,n_sp)
    com2<-rep(0,n_sp)
    names(com1) <- gamma
    names(com2) <- gamma
    ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
    com1[sample(gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1

    ##same for com2:
    com2[sample(gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1

    ##how many species are shared in common?
    #null_shared_spp[i]<-sum((com1+com2)>1)

    ##有两种方法，测试证明，两者结果相似，mantel test结果R为99.1%，p为0.001.这里选择方法2。
    ###方法一：用所有样品平均多度作为概率抽取物种多度（会出现一部分物种多度为0）

    #allreads1 <- sum(spXsite[a1,])
    #com1list <- which(com1==1)
    #com1_pro <- com_sum[com1list]
    #ind1 <- as.factor(sample(com1list,allreads1,replace=TRUE, prob=com1_pro))
    #ind1_summary <- t(as.matrix(table(ind1)))

    #allreads2 <- sum(spXsite[a2,])
    #com2list <- which(com2==1)
    #com2_pro <- com_sum[com2list]
    #ind2 <- as.factor(sample(com2list,allreads2,replace=TRUE, prob=com2_pro))
    #ind2_summary <- t(as.matrix(table(ind2)))


    #pairma <- matrix(0,nrow=2,ncol=n_sp)
    #colnames(pairma) <- 1:n_sp
    #pairma[1,colnames(ind1_summary)] <- ind1_summary
    #pairma[2,colnames(ind2_summary)] <- ind2_summary

    #dis1 <- vegdist(pairma,method="bray")

    ###方法二：以所有样品平均多度为比例，计算物种多度（保证了所有样品多度不为0）
    allreads1 <- sum(spXsite[a1,])
    com1list <- which(com1==1)
    com1_pro <- com_sum[com1list]
    com1_pro <- allreads1*(com1_pro/sum(com1_pro))

    allreads2 <- sum(spXsite[a2,])
    com2list <- which(com2==1)
    com2_pro <- com_sum[com2list]
    com2_pro <- allreads2*(com2_pro/sum(com2_pro))


    pairma <- matrix(0,nrow=2,ncol=n_sp)
    colnames(pairma) <- gamma
    pairma[1,names(com1_pro)] <- com1_pro
    pairma[2,names(com2_pro)] <- com2_pro

    dis2 <- vegdist(pairma,method="bray")


    return(dis2)
  }



  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(apply(spXsite_standard, MARGIN=1, FUN=sum)))

  ##make_null:

  ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels.

  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1

  ##null_array will hold the actual null distribution values.  Each element of the array corresponds to a null distribution for each combination of alpha values.  The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  Later the function will find the row of alpha_table with the right combination of alpha values.  That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels.
  null_array<-list()

  ##looping over each combination of alpha levels:
  
  for(b1 in startn:endn){
    for(b2 in b1:length(alpha_levels)){

      cl <- makeCluster(nocore)
      results <- parLapply(cl,1:reps,func_subRC,spXsite <- spXsite,alpha_levels<-alpha_levels,occur <- occur,a1 <- b1,a2 <- b2)
      stopCluster(cl)
      null_shared_spp <- as.matrix(rbind(results))


      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp

      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[b1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[b2]

      #increment the counter for the columns of the alpha table/ elements of the null array
      print(paste(b1,b2,sep = '_'))

      col_count<-col_count+1

    }
  }

  ##create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")


  #####################
  ##do the test:



  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  numrow <- endn - startn
  results<-matrix(data=NA, nrow=numrow, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))


  ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:numrow){
    for(j in 1:n_sites){

      ##how many species are shared between the two sites:
      #n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      n_shared_obs <- vegdist(spXsite[c(i,j),],method="bray")

      ## what was the observed richness of each site?
      obs_a1<-sum(spXsite_standard[i,])
      obs_a2<-sum(spXsite_standard[j,])

      ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))

      ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))

      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)

      ##how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)



      rc<-(num_greater_in_null)/reps




      if(split_ties){

        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }



      if(!classic_metric){

        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1

        rc<-(rc-.5)*2
      }


      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }

      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }


      ##store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)


    }
  }


  if(as.distance.matrix){
    results<-as.dist(results)
  }


  return(results)

}






