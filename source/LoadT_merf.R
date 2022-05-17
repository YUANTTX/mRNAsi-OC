##
##  标题：把表达数据与临床数据按照相同样本合并；并按照类型聚类
##  作者：  YHJ
##  时间：2021.11.3
##  

clinc_merge <- function(indat,inclinc,intype){
  dat = as.data.frame(indat)
  clinc = as.data.frame(inclinc)
  inset_samples <- intersect(colnames(dat),clinc$submitter)
  dat <- dat[,colnames(dat) %in% inset_samples]
  clinc <- clinc[clinc$submitter %in% inset_samples,]
  
  clinc <- clinc[order(clinc[,intype]),]
  return(clinc)
}

dat_merge <- function(indat,inclinc,intype){
  dat = as.data.frame(indat)
  clinc = as.data.frame(inclinc)
  inset_samples <- intersect(colnames(dat),clinc$submitter)
  dat <- dat[,colnames(dat) %in% inset_samples]
  clinc <- clinc[clinc$submitter %in% inset_samples,]
  
  clinc <- clinc[order(clinc[,intype]),]
  dat <- dat[,clinc$submitter]
  return(dat)
}

