## 
##    step8 建立模型
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


## 1.加载数据
load("rdata/step_1_TCGA_EXP.Rdata")
load("rdata/step_5_TCGA_clinc.Rdata")
load("rdata/step_4_DEG_INFO.Rdata")


## 2.整理数据
DEG_INFO <- DEG_INFO[DEG_INFO$State != 'not',]
COX_dat <- TCGA_EXP
COX_clinc <- TCGA_clinc
COX_dat <- COX_dat[rownames(COX_dat) %in% DEG_INFO$SYMBOL,]

COX_tdat <- as.data.frame(t(COX_dat))
COX_tdat$submitter <- rownames(COX_tdat)
TCGA_clinc <- inner_join(COX_clinc,COX_tdat,by = "submitter")
rownames(TCGA_clinc) <- TCGA_clinc$submitter
colnames(TCGA_clinc) <- gsub(colnames(TCGA_clinc), pattern = '[-]', replacement = '_')
rownames(COX_dat) <- gsub(rownames(COX_dat), pattern = '[-]', replacement = '_')


## 单因素cox
uncox <- uncoxf(indat = TCGA_clinc,ingene = rownames(COX_dat))
uncox <- uncox[uncox$pvalue < 0.05,]
uncox

## lasso回归
lasso <- lassof(indat = TCGA_clinc,ingene = rownames(uncox),insep = 8)
lasso

# 保存数据
write.csv(uncox,file = 'output/step_8_cox_uncox.csv',quote = F)
write.csv(lasso,file = 'output/step_8_cox_lasso.csv',quote = F)
save(lasso,file = 'rdata/step_8_lasso.Rdata')


## 计算风险分数
names(TCGA_clinc)
colnames(TCGA_clinc) <- gsub(colnames(TCGA_clinc), pattern = '-',replacement = '_')
lasso$gene <- gsub(lasso$gene, pattern = '-', replacement = '_')
TCGA_clinc <- theriskf(indat = TCGA_clinc,ingene = lasso$gene,incoeff = lasso$coef)
median(TCGA_clinc$RiskScore)
table(TCGA_clinc$RiskGroup)
save(TCGA_clinc,file = 'rdata/step_8_TCGA_clinc.Rdata')

## ROC绘图
source('source/rocf.R')
rocf(indat = TCGA_clinc)

## 高低风险组KM
sfit <- survfit(Surv(OS.time,OS) ~ RiskGroup, data = TCGA_clinc)
plot4 <- ggsurvplot(sfit,palette = c("#ea0102","#63b931"),
                    risk.table = F,
                    pval =TRUE,pval.size= 4.3,pval.coord=c(1,0.03),
                    pval.method = F,
                    xlim = c(0,150),
                    conf.int = T,xlab ="Time (months)",ylab = 'Survival probability',
                    ggtheme = theme_test(),
                    legend = c(0.82,0.83),
                    legend.title = "Risk Group",legend.labs = c("High risk", "Low risk"))
cox<-coxph(Surv(OS.time,OS) ~ RiskScore,data = TCGA_clinc)
res<-summary(cox)
plot2 <-   plot4$plot +
  annotate("text",x=11,y=0.24,label=paste("HR = ",signif(res$conf.int[1],2))) +
  annotate("text",x=20,y=0.17,label=paste("95%CI = ",round(res$conf.int[3],2),"-",
                                          signif(res$conf.int[4],2)))+
  annotate("text",x=20,y=0.1,label=paste("Cox P = ",signif(res$coefficients[5],3)))
print(plot2)
pdf(file = 'fig/step_8_Mdule_KM.pdf',width = 4,height = 4)
print(plot2)
dev.off()


##  风险因子关联图
source('source/factorf.R')
factorf(indat = TCGA_clinc,inGene = lasso$gene,insep = '8_TCGA')

