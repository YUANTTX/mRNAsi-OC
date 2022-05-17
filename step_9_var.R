## 
##    step9 使用外部数据验证
##


## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(GEOquery))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
source('source/rocf.R')
source('source/oclr.R')
source('source/BaseTabf.R')
source('source/Staticalf.R')
source('source/Id2Symbolf.R')
source('source/factorf.R')


# 加载需要的数据
lasso <- read.csv('output/step_8_cox_lasso.csv',row.names = 1)
load('input/ICGC.Rdata')


# 准备COX数据
COX_dat <- ICGC_DAT
COX_dat <- COX_dat[rownames(COX_dat) %in% lasso$gene,]
COX_clinc <- ICGC_clinc

COX_tdat <- as.data.frame(t(COX_dat))
COX_tdat$icgc_donor_id <- rownames(COX_tdat)
COX_mdat <- inner_join(COX_clinc,COX_tdat,by = "icgc_donor_id")
rownames(COX_mdat) <- COX_mdat$icgc_donor_id
COX_mdat$OS.time <- as.numeric(COX_mdat$OS.time)
COX_mdat$OS <- as.numeric(COX_mdat$OS)


## 计算分线得分
Risk_mdat <- theriskf(indat = COX_mdat,ingene = lasso$gene,incoeff = lasso$coef)
Risk_mdat$RiskScore
median(Risk_mdat$RiskScore)
table(Risk_mdat$RiskGroup)

## ROC绘图
rocf(indat = Risk_mdat)

## 高低风险组KM
sfit <- survfit(Surv(OS.time,OS) ~ RiskGroup, data = Risk_mdat)
plot4 <- ggsurvplot(sfit,palette = c("#ea0102","#63b931"),
                    risk.table = F,
                    pval =TRUE,pval.size= 4.3,pval.coord=c(1,0.03),
                    pval.method = F,
                    xlim = c(0,150),
                    conf.int = T,xlab ="Time (months)",ylab = 'Survival probability',
                    ggtheme = theme_test(),
                    legend = c(0.82,0.83),
                    legend.title = "Risk Group",legend.labs = c("High risk", "Low risk"))
cox<-coxph(Surv(OS.time,OS) ~ RiskScore,data = Risk_mdat)
res<-summary(cox)
plot2 <-   plot4$plot +
  annotate("text",x=11,y=0.24,label=paste("HR = ",signif(res$conf.int[1],2))) +
  annotate("text",x=20,y=0.17,label=paste("95%CI = ",round(res$conf.int[3],2),"-",
                                          signif(res$conf.int[4],2)))+
  annotate("text",x=20,y=0.1,label=paste("Cox P = ",signif(res$coefficients[5],3)))
print(plot2)
pdf(file = 'fig/step_9_var_KM.pdf',width = 4,height = 4)
print(plot2)
dev.off()

