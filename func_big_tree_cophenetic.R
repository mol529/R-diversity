##############################################################
####以下代码需要在服务器跑，个人电脑可能带不动！！！！########
##############################################################

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



