##
##  标题:诺莫图
##  作者：yhj
##  时间：2021.10.18
## 

suppressMessages(library(rms))
suppressMessages(library(survival))


nomodulef <- function(indat,invar,instep,inmm){

  Nmdat <- indat
  Nmdat <- Nmdat[,c(invar,'OS','OS.time')]
  Nmdat$OS.time <- as.numeric(Nmdat$OS.time)*30
  dd <- datadist(Nmdat)
  options(datadist="dd")
  multivarl <- as.formula(paste0('Surv(OS.time,OS)~', 
                                 paste(invar, sep = '', collapse = '+')))
  
  print('Building the nomogram graph')
  coxm_1 <- cph(formula = multivarl,data=Nmdat,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  nomo <- nomogram(coxm_1,fun = list(surv1,surv3,surv5),lp = T,
                funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
                maxscale = 100,fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  pdf(file = paste0('fig/step_',instep,'_nomodule.pdf'),width = 4,height = 4)
  plot(nomo,lplabel = 'Linear Preadictor',
       xfrac = .35,varname.label = T,varname.label.sep = '=',ia.space = .2,
       tck = NA,tcl = 0.2,lmgp = 0.3,
       points.label = 'Points',total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
  dev.off()
  print('successful')
  print('Analyzing forecast performance')
  f1 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=inmm, B=1000)
  f3 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=inmm, B=1000)
  f5 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=inmm, B=1000)
  
  pdf(file = paste0('fig/step_',instep,'_nomoadjust.pdf'),width = 4,height = 4)
  plot(cal1,lwd=1,lty=1, cex.axis = 1,cex.lab=1,
       errbar.col = '#666699',
       xlab='Nomogram-Predicted Probability',
       ylab='Actual',
       col = '#666699',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal3,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#339933',
       col = '#339933',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal5,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#FF0033',
       col = '#FF0033',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  abline(0,1,lty=1,lwd=1)
  legend("bottomright",legend=c("1 - year","3 - year","5 - year"), 
         col=c("#666699","#339933","#FF0033"),
         # x.intersp=3, y.intersp=0.8,
         lty= 1 ,lwd= 4,
         bty = "n",
         seg.len=1,cex=1)
  dev.off()
  
  c_in <- coxph(multivarl,data=Nmdat)
  sum.surv <- summary(c_in)
  print(sum.surv$concordance)
  print('successful')
  return(nomo)
}

