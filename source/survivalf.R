##
##  标题：生存分析函数
##  作者：YHJ
##  时间：2021.11.3
##  备注：需要提前将基因合并在临床数据中；简单+复杂
##

suppressMessages(library(survival))
suppressMessages(library(survminer))


Surf <- function(indat,invar,inp){
  surdat <- indat
  genelist <- invar
  surTab <- data.frame()
  pFilter <- inp
  dir.create('fig/sur',recursive = T)
  
  for(gene in genelist){
    if(sd(surdat[,gene])==0){
      next}
    a      <- surdat[,gene]<=median(surdat[,gene])
    diff   <- survdiff(Surv(OS.time, OS) ~ a,data = surdat)
    pValue <- 1-pchisq(diff$chisq,df=1)
    surTab <- rbind(surTab,cbind(gene=gene,pvalue=pValue))
    
    SurDiff <- survSurDiff(Surv(OS.time, OS) ~ a, data = surdat)
    if(pValue < pFilter){
      if(pValue < 0.05){
        pValue <- signif(pValue,4)
        pValue <- format(pValue, scientific = TRUE)
      }else{
        pValue <- round(pValue,3)
      }
      pdf(file=paste0("fig/sur/",gene,".survival.pdf"),width = 4,height =4)
      plot(SurDiff, 
           lwd=1,cex = 1,
           cex.axis=1,cex.lab=1,
           col=c("red","blue"),
           xlab="OS.time (mouth)",
           mark.time=T,
           ylab="Survival rate",
           ylim=c(0,1.09),
           main=paste(gene,": Survival curve (p=", pValue ,")",sep=""))
      legend("topright", 
             c(paste(gene," high expression",sep=""), 
               paste(gene," low expression",sep="") ), 
             lwd=0.5, cex=0.5,
             col=c("red","blue"))
      dev.off()
    }
  }
  
  return(surTab)
}

