#
#     标题：使用forestplot绘制本地cox数据
#     来源：本地数据COX数据
#     时间：2021.10.18
#     备注：格式已准备好
#


Forstf <- function(indat,instep){
  suppressMessages(library(forestplot))
  coxdat <- indat
  hrtext <- cbind(c(' ',rownames(coxdat)),
                  c('P.value',coxdat$pvalue),
                  c('Hazard ratio',coxdat$HR_95_CI))
  hrdat <- data.frame('mean'  = c(NA,coxdat$HR),
                      'lower' = c(NA,coxdat$HR95L),
                      'upper' = c(NA,coxdat$HR95H))
  hrdat <- apply(as.matrix(hrdat),2,function(x) as.numeric(x))
  plot1 <-   forestplot(hrtext, hrdat,
                        boxsize = 0.2,
                        clip=c(0,2.5),
                        zero = 1,lwd.zero = 3,lwd.ci = 3,
                        xlab="Hazard ratio",
                        colgap=unit(4,"mm"),
                        txt_gp=fpTxtGp(label = gpar(cex = 1),
                                       ticks = gpar(cex = 1),
                                       xlab  = gpar(cex = 1),
                                       title = gpar(cex = 1)),
                        col=fpColors(box="#00c853",line="darkblue", summary="royalblue"))
  pdf(file = paste0('fig/step_',instep,'_forst.pdf'),width = 5,height = signif(length(rownames(coxdat))/2,2))
  print(plot1)
  dev.off()
}
