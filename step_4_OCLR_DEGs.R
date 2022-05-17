##
##    step4 差异分析及通路富集
## 


## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(dplyr))
source('source/DEGf.R')
source('source/DEGplotf.R')
source('source/pathwayf.R')
source('source/BaseTabf.R')


## 加载数据
load('rdata/step_1_TCGA_EXP.Rdata')
load('rdata/step_1_TCGA_clinc.Rdata')


## mRNAsi生存分析确定最佳的截断值,并添加到TCGA_clinc
names(TCGA_clinc)
TCGA_clinc[,c(100,101)] <- char2numf(TCGA_clinc[,c(100,101)])
rownames(TCGA_clinc) <- TCGA_clinc$submitter

res.cut <- surv_cutpoint(TCGA_clinc,
                         time = "OS.time", 
                         event = "OS", 
                         variables = c("mRNAsi"))
plot(res.cut, "mRNAsi", palette = "npg")
pdf(file = 'fig/step_4_mRNAsi_group.pdf',width = 4,height = 5)
plot(res.cut, "mRNAsi", palette = "npg")
dev.off()

res.cat <- surv_categorize(res.cut)
res.cat$submitter <- rownames(res.cat)
res.cat$mRNAsiGroup2 <- res.cat$mRNAsi
res.cat <- res.cat[,c('submitter','mRNAsiGroup2')]
TCGA_clinc <- inner_join(TCGA_clinc,res.cat, by = 'submitter')

diff = survdiff(Surv(OS.time, OS) ~ mRNAsiGroup2,data = TCGA_clinc)
pValue = 1 - pchisq(diff$chisq, df = 1)
sfit <- survfit(Surv(OS.time,OS) ~ mRNAsiGroup2, data = TCGA_clinc)
plot4 <- ggsurvplot(sfit,palette = c("#248888","#ea7070"),
                    risk.table = T,
                    pval =TRUE,pval.size= 4.3,pval.coord=c(1,0.03),
                    pval.method = F,
                    conf.int = F,xlab ="Time (months)",ylab = 'Survival probability',
                    ggtheme = theme(),
                    legend = c(0.82,0.83),
                    legend.title = "",legend.labs = c( "Low mRNAsi","High mRNAsi"))
cox<-coxph(Surv(OS.time,OS) ~ mRNAsiGroup2,data = TCGA_clinc)
res<-summary(cox)
plot2 <-   plot4$plot +
  annotate("text",x=24,y=0.29,label=paste("HR = ",signif(res$conf.int[1],2))) +
  annotate("text",x=38,y=0.22,label=paste("95%CI = ",round(res$conf.int[3],2),"-",
                                          signif(res$conf.int[4],2)))+
  annotate("text",x=36,y=0.15,label=paste("Cox P = ",signif(res$coefficients[5],3)))
print(plot2)
pdf(file = 'fig/step_4_mRNAsi_KM2.pdf',width = 4,height = 5)
print(plot2)
dev.off()


## 样本分组
TCGA_clinc <- TCGA_clinc[order(TCGA_clinc$mRNAsiGroup2),]
group_mRNAsi <- TCGA_clinc$mRNAsiGroup2
TCGA_EXP <- TCGA_EXP[,TCGA_clinc$submitter]


## 差异分析
batoh = 'high-low'
DEG_INFO <- limmaf(indat = TCGA_EXP,ingroup = group_mRNAsi,
                   batoh = batoh,inp = 0.05,inlogfc = 1)
table(DEG_INFO$State)


## 差异结果绘图
DEGPlotf(indat = TCGA_EXP,ininfo = DEG_INFO,intype = c('high','low'),
         ingroup = group_mRNAsi,insep = 4)


## 通路富集分析
Pay_INFO <- DEG_INFO[DEG_INFO$State != 'not',]
Pathwayf(ingene = Pay_INFO$SYMBOL,insep = 4)


## 保存数据
save(TCGA_clinc,file = 'rdata/step_4_TCGA_clinc.Rdata')
save(DEG_INFO,file = 'rdata/step_4_DEG_INFO.Rdata')

