##
##  标题：1.计算分线得分 2.绘制1、3、5年ROC曲线
##  作者：YHJ
##  时间：2021.11.3
##  

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(survivalROC))


theriskf <- function(indat,ingene,incoeff){
  
  ScoreDat <- indat
  filename <- ingene
  beta <- incoeff
  
  print('Calculating score')
  dat_score <- 0    
  for (i in 1:length(filename)){
    fm <- filename[i]
    bt <- beta[i]
    bt <- as.numeric(bt)
    sc <- ScoreDat[,fm]*bt
    dat_score <- dat_score + sc
  }
  ScoreDat$RiskScore <- dat_score
  ScoreDat$RiskGroup <- ifelse(ScoreDat$RiskScore > median(ScoreDat$RiskScore,na.rm = T),'high','low')
  print('successful')
  
  return(ScoreDat)
}

rocf <- function(indat){
  
  ROCDat <- indat
  print('ROC analysis in progress')
  cuto <- cbind(12, 36, 60, 84)
  nobs <- NROW(ROCDat)
  for (i in 1:length(cuto)){
    cutoff <- cuto[i]
    # Mayo4 <- paste('Mayo4',cuto[i],sep = '.')
    Mayo4= survivalROC(Stime=ROCDat$OS.time,  
                       status=ROCDat$OS,
                       marker = ROCDat$RiskScore,     
                       predict.time =  cutoff, 
                       method="KM",
                       span = 0.25*nobs^(-0.20))
    plot(Mayo4$FP, Mayo4$TP, lwd=2,
         type="l",col="#FF6666",xlim=c(0,1), ylim=c(0,1),   
         xlab='False positive rate',
         ylab="True positive rate")
    abline(0,1,col="black",lwd=1)
    cut.op= Mayo4$cut.values[which.max(Mayo4$TP-Mayo4$FP)]
    print(cut.op)
    Mayo4$AUC
    assign(paste0('Mayo4','.',cutoff),Mayo4) 
  }
  
  lines(Mayo4.12$FP, Mayo4.12$TP, type="l",col="#FFFF66",xlim=c(0,1), ylim=c(0,1),lwd=2)
  lines(Mayo4.36$FP, Mayo4.36$TP, type="l",col="#66CC00",xlim=c(0,1), ylim=c(0,1),lwd=2)
  lines(Mayo4.60$FP, Mayo4.60$TP, type="l",col="#0099CC",xlim=c(0,1), ylim=c(0,1),lwd=2)
  legend(0.45,0.3,c(paste("AUC at 1 years:",round(Mayo4.12$AUC,3)),
                    paste("AUC at 3 years:",round(Mayo4.36$AUC,3)),
                    paste("AUC at 5 years:",round(Mayo4.60$AUC,3)),
                    paste("AUC at 7 years:",round(Mayo4.84$AUC,3))),
         x.intersp=1, y.intersp=0.7,
         lty= 1 ,lwd= 2,col=c('#FFFF66',"#66CC00","#0099CC",'#FF6666'),
         bty = "n",
         seg.len=1,cex=0.65)
  print('sucessful')
}

# 批量绘制ROC
# library(ggsci)
# library(survivalROC)
# color2 = c('#F97F51','#B33771','#6D214F','#EAB543','#1B9CFC','#FD7272',
#            '#FC427B','#55E6C1','#58B19F','#BDC581','#D6A2E8','#F97F51')
# ROC_mclinc <- TCGA_clinc
# varlist = list()
# 
# pdf(file = 'fig/step_9_genes_ROC.pdf',width = 7,height = 7)
# plot(0, 0, type = "l",xlim = c(0,1), ylim = c(0,1),
#      xlab = "False positive rate", ylab = "True positive rate")
# abline(0,1,col="black")
# for (i in 1:length(lasso$gene)) {
#   k = lasso$gene[i]
#   nobs <- NROW(ROC_mclinc)
#   Mayo4 = survivalROC(Stime=ROC_mclinc$OS.time,  
#                       status=ROC_mclinc$OS,      
#                       marker = ROC_mclinc[,k],
#                       predict.time =  72,
#                       span = 0.25*nobs^(-0.20))
#   assign('Mayo4_tmp',Mayo4) 
#   lines(Mayo4_tmp$FP, Mayo4_tmp$TP, type="l",col=color2[i],lwd=1.5)
#   
#   tmplist = list(paste("AUC of ",k," =",round(Mayo4_tmp$AUC,3)))
#   varlist = c(tmplist,varlist)
# }
# legend("bottomright",unlist(varlist),
#        x.intersp = 1, y.intersp = 0.9,
#        lty = 1 ,lwd = 1.5,col = color2,
#        bty = "o",
#        seg.len = 1,cex = 0.75)
# dev.off()