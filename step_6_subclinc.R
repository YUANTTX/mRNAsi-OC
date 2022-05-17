##
##    step6 干性新亚型评估：
##


## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggsci))
suppressMessages(library(stringr))
suppressMessages(library(maftools))
suppressMessages(library(dplyr))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(limma))
source('source/BaseTabf.R')


## 加载数据
load('rdata/step_1_TCGA_EXP.Rdata')
load('rdata/step_5_TCGA_clinc.Rdata')
load('rdata/step_4_DEG_INFO.Rdata')

TCGA_clinc <- TCGA_clinc[order(TCGA_clinc$cluster),]
TCGA_EXP <- TCGA_EXP[,TCGA_clinc$submitter]
TCGA_EXP[,1:ncol(TCGA_EXP)] <- char2numf(TCGA_EXP[,1:ncol(TCGA_EXP)])

## 绘制热图
HtDat <- TCGA_EXP[DEG_INFO[DEG_INFO$State != 'not',]$SYMBOL,]
annotation_col = data.frame(mRNAsi = TCGA_clinc$mRNAsi,
                            ImmuneScore = TCGA_clinc$ImmuneScore,
                            Steamness_Subtype = TCGA_clinc$StemType)
rownames(annotation_col) = TCGA_clinc$submitter
ann_colors = list(Steamness_Subtype = c('Steamness Subtype I'="#e19d29", 'Steamness Subtype II'="#4298b5"),
                  ImmuneScore = c('#b3e5fc','#0277bd'),
                  ESTIMATEScore = c('#fff9c4','#f57f17'),
                  TumorPurity = c('#c8e6c9','#2e7d32'),
                  mRNAsi = c('#e1bee7','#6a1b9a'))
mian_col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
plot8 <- pheatmap(as.matrix(HtDat), name = " ",
                  cluster_cols = F,cluster_rows = T,
                  annotation_col = annotation_col, 
                  annotation_colors = ann_colors,
                  treeheight_row = 8,
                  show_colnames = F,show_rownames = F,
                  scale="row",col = mian_col)
plot8
pdf(file="fig/step_5_Stemtype_heatmap.pdf",width = 8,height = 6)
plot8
dev.off()


## 亚组生存分析
names(TCGA_clinc)
TCGA_clinc[,c(100,101)] <- char2numf(TCGA_clinc[,c(100,101)])

sfit <- survfit(Surv(OS.time,OS) ~ StemType, data = TCGA_clinc)
plot4 <- ggsurvplot(sfit,palette = c("#e19d29","#4298b5"),
                    risk.table = F,
                    pval =TRUE,pval.size= 4.3,pval.coord=c(1,0.03),
                    pval.method = F,
                    xlim = c(0,150),
                    conf.int = T,xlab ="Time (months)",ylab = 'Survival probability',
                    ggtheme = theme_test(),
                    legend = c(0.82,0.83),
                    legend.title = "Subtype",
                    legend.labs = c("Subtype I", "Subtype II"))
pdf(file = 'fig/step_6_Stemtype_km.pdf',width = 5,height = 5)
print(plot4)
dev.off()

## 亚组GSVA分析
geneset <-  getGmt('input/h.all.v7.0.symbols.gmt')
GSVADat <- gsva(as.matrix(TCGA_EXP), geneset, method = "gsva",ssgsea.norm=TRUE,verbose=TRUE)

design <- model.matrix(~ factor(TCGA_clinc$StemType))
colnames(design) <- c("Steamness Subtype I", "Steamness Subtype II")
row.names(design)<-colnames(TCGA_EXP)
fit <- lmFit(GSVADat, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="Steamness Subtype II", number=Inf)
DEgeneSets <- topTable(fit, coef="Steamness Subtype II", number=Inf,
                       p.value = 0.05, adjust="BH")

GSVADat <- GSVADat[rownames(GSVADat) %in% rownames(DEgeneSets),]
rownames(GSVADat) <- gsub(rownames(GSVADat) , pattern = 'HALLMARK_', replacement = '')
annotation_col = data.frame(StemType = TCGA_clinc$StemType)
rownames(annotation_col) = TCGA_clinc$submitter
ann_colors = list(StemType = c('Steamness Subtype I'='#e19d29','Steamness Subtype II'='#4298b5'))
mian_col = colorRamp2(c(-2, 0, 2), c("blue", "#263238", "#ffd600"))
plot1 <- pheatmap(as.matrix(GSVADat), name = " ",
                  cluster_cols = F,cluster_rows = T,
                  annotation_col = annotation_col, 
                  annotation_colors = ann_colors,
                  treeheight_row = 8,annotation_legend = F,
                  show_colnames = F,show_rownames = T,
                  scale="row",col = mian_col,
                  border_color = F)
pdf(file="fig/step_6_StemType_GSVA.pdf",width = 8,height = 6)
plot1
dev.off()


##  突变负荷TMB比较
lamlall <- read.maf(maf = "input/OC_Maf.txt")
Maf_Tmb = tmb(maf = lamlall, captureSize = 50, logScale = T)
Maf_Tmb$submitter <- substr(Maf_Tmb$Tumor_Sample_Barcode,1,16)
Maf_Tmb$submitter <- gsub(Maf_Tmb$submitter, pattern = '[-]', replacement = '_')
TCGA_clinc2 <- inner_join(TCGA_clinc,Maf_Tmb,by = 'submitter')
color2 = c('#e19d29','#4298b5')
comparison = list(c("Steamness Subtype I","Steamness Subtype II"))
names(TCGA_clinc2)
boxp9 = ggboxplot(TCGA_clinc2, size = 1, bxp.errorbar = F,
                  x="StemType", y="total_perMB", color = "StemType",
                  palette = color2, ylab = 'TMB',
                  order = c("Steamness Subtype I","Steamness Subtype II")) + 
  stat_compare_means(comparisons = comparison, method = "wilcox.test") +
  theme(text = element_text(size = 10),legend.text = element_text(size = 14),
        legend.position = 'none') + 
  scale_x_discrete('StemType',labels = c("Steamness Subtype I" = "Subtype I","Steamness Subtype II" = "Subtype II"))
boxp9
ggsave(plot = boxp9,filename = "fig/step_6_Stemtype_TMB.pdf",width = 3,height = 4)


