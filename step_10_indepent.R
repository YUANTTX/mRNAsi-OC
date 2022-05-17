## 
##    step10 预后模型独立性 + nomogram
## 

## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
source('source/rocf.R')
source('source/BaseTabf.R')
source('source/Staticalf.R')


## 加载cox结果
load('rdata/step_8_lasso.Rdata')
load('rdata/step_8_TCGA_clinc.Rdata')

## 整理数据
names(TCGA_clinc)
colnames(TCGA_clinc) <- gsub(colnames(TCGA_clinc), pattern = '-',replacement = '_')
lasso$gene <- gsub(lasso$gene, pattern = '-', replacement = '_')

## 计算得分
TCGA_clinc <- theriskf(indat = TCGA_clinc,ingene = lasso$gene,incoeff = lasso$coef)
median(TCGA_clinc$RiskScore)

## 临床相关性分析
names(TCGA_clinc)
invarlist <- c("RiskScore",'Stage2','Age2')
TCGA_clinc$Age2

uncox <- uncoxf(indat = TCGA_clinc,ingene = invarlist)
uncox
mucox <- mucoxf(indat = TCGA_clinc,ingene = invarlist)
mucox

## 绘制森林图
library(forestplot)
coxlist <- c('uncox','mucox')
collist <- c('green','red')
for (i in 1:2) {
  indat <- get(coxlist[i])
  hrtext <- cbind(c(' ',rownames(indat)),
                  c('P.value',indat$pvalue),
                  c('Hazard ratio',indat$HR_95_CI))
  hrdat <- data.frame('mean'  = c(NA,indat$HR),
                      'lower' = c(NA,indat$HR95L),
                      'upper' = c(NA,indat$HR95H))
  hrdat <- apply(as.matrix(hrdat),2,function(x) as.numeric(x))
  plot1 <- forestplot(hrtext, hrdat,
                      boxsize = 0.2,
                      zero = 1,lwd.zero = 1,lwd.ci = 1,
                      xlab="Hazard ratio",
                      graphwidth = unit(4,'cm'),
                      lineheight = unit(1.5,'cm'),
                      colgap=unit(4,"mm"),
                      txt_gp=fpTxtGp(label = gpar(cex = 1),
                                     ticks = gpar(cex = 1),
                                     xlab  = gpar(cex = 1),
                                     title = gpar(cex = 1)),
                      col=fpColors(box=collist[i],line="darkblue", summary="royalblue")
  )
  print(plot1)
  pdf(file = paste0('fig/step_10_Independ_',coxlist[i],'.pdf'),width = 6,height = 4)
  print(plot1)
  dev.off()
}


## 模型基因的表达热图
load('rdata/step_4_DEG_INFO.Rdata')

TCGA_clinc <- TCGA_clinc[order(TCGA_clinc$RiskScore),]
htdat <- TCGA_clinc[,c(lasso$gene)]
htdat <- as.data.frame(t(htdat))
rownames(DEG_INFO) <- gsub(rownames(DEG_INFO) , pattern = '-', replacement = '_')
DEG_INFO2 <- DEG_INFO[rownames(DEG_INFO) %in% rownames(htdat),]

x = DEG_INFO2$logFC
names(x) = rownames(DEG_INFO2)
FCgene = names(sort(x))
ht_dat = htdat[FCgene,]

ann_col = data.frame(StemType   = TCGA_clinc$StemType,
                     Age        = TCGA_clinc$Age2,
                     Stage      = TCGA_clinc$Stage2,
                     RiskGroup  = TCGA_clinc$RiskGroup)
rownames(ann_col) = colnames(ht_dat)
ann_colors = list(Age       = c('>60'="#FF9999",'<=60'="#FFCC99"),
                  Stage     = c("I/II"="#99CCCC","III/IV"="#FFCC99"),
                  StemType  = c('Steamness Subtype I' = '#e19d29','Steamness Subtype II' = '#4298b5'),
                  RiskGroup = c('low' = '#63b931','high' = '#ea0102'))
main_col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
plot6 <- pheatmap(as.matrix(ht_dat),name = ' ',
                  show_rownames = T,show_colnames = F,
                  cluster_rows=F,cluster_cols = F,
                  annotation_row = NA,
                  annotation_col = ann_col,
                  annotation_colors =ann_colors,
                  treeheight_row = 8,
                  scale="row",color = main_col)
plot6
pdf(file="fig/step_10_Independ_heatmap.pdf",width = 6,height = 4)
plot6
dev.off()


## 诺莫图
indat = TCGA_clinc
invar = c('RiskScore')
instep = 10
inmm = 100
{
  library(rms)
  library(survival)
  Nmdat <- indat
  Nmdat <- Nmdat[,c(invar,'OS','OS.time')]
  Nmdat$OS.time <- as.numeric(Nmdat$OS.time)*30
  dd <- datadist(Nmdat)
  options(datadist="dd")
  multivarl <- as.formula(paste0('Surv(OS.time,OS)~', 
                                 paste(invar, sep = '', collapse = '+')))
  
  coxm_1 <- cph(formula = multivarl,data=Nmdat,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  nomo <- nomogram(coxm_1,fun = list(surv1,surv3,surv5),lp = T,
                   funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
                   maxscale = 100,fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  pdf(file = paste0('fig/step_',instep,'_Module_nomodule.pdf'),width = 4,height = 4)
  plot(nomo,lplabel = 'Linear Preadictor',
       xfrac = .35,varname.label = T,varname.label.sep = '=',ia.space = .2,
       tck = NA,tcl = 0.2,lmgp = 0.3,
       points.label = 'Points',total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
  dev.off()
  ## 校准曲线
  f1 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=inmm, B=1000)
  f3 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=inmm, B=1000)
  f5 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=inmm, B=1000)
  
  pdf(file = paste0('fig/step_',instep,'_Module_nomoadjust.pdf'),width = 4,height = 4)
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
         lty= 1 ,lwd= 4,
         bty = "n",
         seg.len=1,cex=1)
  dev.off()
}
