library(picante)
source("func_big_tree_cophenetic.R")
source("func_bNTI.R")
##读取otu table
featuretable <- "rarefied-9622-table.txt"
com <- read.table(featuretable,header=T,row.names = 1, sep='\t',check.name = F)
com <- t(com)

grp <- "group_yr.txt"
group <- read.table(grp,header=T,sep='\t',check.name = F)

##读取系统发育�?
tree<-read.tree(file = "rooted-tree.nwk")
tree1 <- prune.sample(com,tree)
cop_all <- bigtree_cophenetic(tree1)



##分类计算
for (n in c('AWP','AWF','AS','SWP','SWF','SS')) {
  sam_list <- as.character(group[group$type==n,"sampleid"])
  com_sub <- com[sam_list,]
  spe_list <- colnames(com_sub)[which(colSums(com_sub)>0)]
  
  com_sub <- com_sub[,spe_list]
  cop_sub <- cop_all[spe_list,spe_list]
  betaNTI(com_sub,cop_sub,filename = n,nocore = 24,beta.reps = 999)

}


##所有样品一起计�?
betaNTI(com,cop_all,filename = "All",nocore = 24,beta.reps = 999)
