##
##  标题：一致性聚类:只需要一个表达矩阵
##  作者：YHJ
##  时间：2021.11.3
##

suppressMessages(library(ConsensusClusterPlus))

newcluster <- function(indat,ink){
  
  print('Unsupervised clustering in progress')
  if (!file.exists('cluster/fig')){
    dir.create('cluster/fig',recursive = T)
  }
  sub_dat <- as.matrix(indat)
  mads <- apply(sub_dat,1,mad)
  sub_dat <-  sweep(sub_dat,1, apply(sub_dat,1,median,na.rm=T))
  maxK <-  ink
  results <-  ConsensusClusterPlus(sub_dat,
                                   maxK = maxK,
                                   reps = 50,
                                   seed = 1,
                                   pItem = 0.8,
                                   pFeature = 1,
                                   clusterAlg = "pam",
                                   distance="pearson",
                                   title="./cluster",
                                   innerLinkage="complete",
                                   writeTable = T,
                                   plot="png")
  Kvec = 2:maxK
  x1 = 0.1; x2 = 0.9
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="")
  for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  } 
  optK = Kvec[which.min(PAC)]
  print(paste0('The best ink is: ',optK))
  icl = calcICL(results,title="./cluster/fig/",plot="pdf")
  print('successful')
}