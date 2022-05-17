## 
##    step1: 数据整理：表达数据，临床数据统一格式
##    时间：2021.8.26
##    备注：
##    肿瘤：卵巢癌
##
 
## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(pacman))
source('source/LoadT_Datf.R')
source('source/oclr.R')
source('source/immunef.R')
source('source/fpkm2Tpm.R')


## 加载数据: 
TCGA_EXP <- Readdat(inpath = 'input/TCGA-OV.htseq_fpkm.tsv.gz',
                     idpath = 'input/gencode.v22.annotation.gene.probeMap')
TCGA_clinc <- ReadClinc(inpath = 'input/TCGA-OV.GDC_phenotype.tsv.gz')

#分别保存OCLR与ESTIMATE的读入文件
TCGA_oclr <- cbind(gene_id = rownames(TCGA_EXP),TCGA_EXP)
write.table(TCGA_oclr,file = 'tmp/step_1_TCGA_oclr.txt',sep = '\t',quote = F)
write.table(TCGA_EXP,file = 'tmp/step_1_TCGA_EXP.txt',sep = '\t',quote = F)


## 临床数据:提取有用的临床信息
TCGA_clinc$OS.time <- as.numeric(signif(TCGA_clinc$days_to_last_follow_up.diagnoses/30,3))
TCGA_clinc$OS <- ifelse(TCGA_clinc$vital_status.demographic == 'Alive',0,
                        ifelse(TCGA_clinc$vital_status.demographic == 'Dead',1,'NA'))
TCGA_clinc$OS <- as.numeric(TCGA_clinc$OS)
TCGA_clinc$Status <- ifelse(TCGA_clinc$OS == 0,'Alive',ifelse(TCGA_clinc$OS == 1,'Dead','NA'))
TCGA_clinc$Age <- TCGA_clinc$age_at_initial_pathologic_diagnosis
TCGA_clinc$Age2 <- ifelse(TCGA_clinc$Age > 60,'>60','<=60')
TCGA_clinc$Stage <- str_extract(TCGA_clinc$clinical_stage,pattern = "[IV]+")
TCGA_clinc$Stage2 <- ifelse(TCGA_clinc$Stage %in% c('I','II'),'I/II','III/IV')
TCGA_clinc$histologic_grade <- TCGA_clinc$neoplasm_histologic_grade
TCGA_clinc$residual <- TCGA_clinc$tumor_residual_disease


## 计算干性分数：StenScore
mRNAsi <- oclrpd(sigpath = 'input/pcbc-stemsig.tsv',
                  datpath = 'tmp/step_1_TCGA_oclr.txt')

## 计算免疫分数: ESTIMATE
ESTIscore <- Immestif(datpath = 'tmp/step_1_TCGA_EXP.txt',inplat = 'illumina',insep = 1)

## 整理数据
TCGA_clinc <- Reduce(function(x,y) inner_join(x,y,by="submitter",all.x=TRUE),
                     list(TCGA_clinc,mRNAsi,ESTIscore),accumulate =FALSE)

## 保存数据
save(TCGA_EXP, file = 'rdata/step_1_TCGA_EXP.Rdata')
save(TCGA_clinc,file = 'rdata/step_1_TCGA_clinc.Rdata')

