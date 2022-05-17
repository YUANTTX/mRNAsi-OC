## 
##    step5 根据mRNAsi差异分析结果分亚群
## 

## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(dplyr))
source('source/subclusterf.R')

## 加载数据
load('rdata/step_1_TCGA_EXP.Rdata')
load('rdata/step_4_TCGA_clinc.Rdata')
load('rdata/step_4_DEG_INFO.Rdata')

##  整理数据
DEG_dat <- TCGA_EXP[DEG_INFO[DEG_INFO$State != 'not',]$SYMBOL,]

##  无监督聚类
newcluster(indat = DEG_dat,ink = 10)

##  添加到临床信息
Kcluster <- read.csv('cluster/cluster.k=2.consensusClass.csv',header = F)
names(Kcluster) <- c('submitter','cluster')

TCGA_clinc <- left_join(TCGA_clinc,Kcluster,by = 'submitter')
TCGA_clinc$StemType <- ifelse(TCGA_clinc$cluster == 1,'Steamness Subtype I','Steamness Subtype II')

##  保存数据
save(TCGA_clinc,file = 'rdata/step_5_TCGA_clinc.Rdata')




