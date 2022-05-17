##
##    step3 干性指数与免疫关系
## 


## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(GSVA))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(ggpubr))
suppressMessages(library(ggExtra))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(rlang))
suppressMessages(library(ggside))
source('source/immunef.R')
source('ImmDAT/CIBERSORT/Cibersort.R')

## 加载数据
load('rdata/step_1_TCGA_EXP.Rdata')
load('rdata/step_1_TCGA_clinc.Rdata')
load("ImmDAT/imm1316/ImmGene_Set.Rdata")


##  计算ssGSEA免疫分数
ssGSEAres <- ImmssGSEAf(indat =TCGA_EXP,inSet = gene.set.list)
TCGA_clinc <- TCGA_clinc[order(TCGA_clinc$mRNAsi),]
ssGSEAres <- ssGSEAres[,TCGA_clinc$submitter]

## 绘制ssGSEA免疫热图
annotation_col = data.frame(StromalScore  = TCGA_clinc$StromalScore,
                            ImmuneScore   = TCGA_clinc$ImmuneScore,
                            ESTIMATEScore = TCGA_clinc$ESTIMATEScore,
                            mRNAsiGroup     = TCGA_clinc$mRNAsiGroup)
rownames(annotation_col) = TCGA_clinc$submitter
ann_colors = list(StromalScore  = c("white", "firebrick"),
                  ImmuneScore   = c('#e0f7fa','#26c6da'),
                  ESTIMATEScore = c('#fff9c4','#f57f17'),
                  mRNAsiGroup     = c('low'='#248888','high'='#ea7070'))
mian_col = colorRamp2(c(-4, 0, 4), c("#515bd4", "white", "red"))
plot1 <- pheatmap(as.matrix(ssGSEAres), name = "Immune",
                  cluster_cols = F,cluster_rows = T,
                  annotation_col = annotation_col, 
                  annotation_colors = ann_colors,
                  show_colnames = F,show_rownames = T,
                  scale="row",col = mian_col)
plot1
pdf(file="fig/step_3_mRNAsi_ssGSEA.pdf",width = 12,height = 6)
plot1
dev.off()


## OCLR与免疫分数相关性
score = c('StromalScore','ImmuneScore','ESTIMATEScore')
for (i in score) {
  invar <- sym(i)
  plot2 <- ggscatterstats(TCGA_clinc, 
                          x = !!invar, y = "mRNAsi",
                          bf.message = F,marginal.type = 'density',
                          ggtheme = theme_bw())
  print(plot2)
  pdf(file = paste0('fig/step_3_mRNAsi_',i,'.pdf'),width = 4,height = 4)
  print(plot2)
  dev.off()
}

## OCLR与22种免疫细胞相关性
LM22.file <- "ImmDAT/CIBERSORT/LM22.txt"
TCGA_exp.file = 'tmp/step_1_TCGA_EXP.txt'

TCGA_TME.results <- CIBERSORT(LM22.file, TCGA_exp.file, perm = 1000, QN = F)
TCGA_TME.results <- TCGA_TME.results[TCGA_clinc$submitter,]

group_list <- TCGA_clinc$mRNAsiGroup %>% factor(.,levels = c("low","high"))

TME_data <- as.data.frame(TCGA_TME.results[,1:22])
TME_data$group <- group_list
TME_data$sample <- row.names(TME_data)
TME_New = melt(TME_data)
colnames(TME_New)=c("Group","Sample","Celltype","Immune cell composition")

boxpplot <- ggboxplot(TME_New, x = "Celltype", y = "Immune cell composition",
                      color = "Group", palette = c('#248888','#ea7070')) + 
  stat_compare_means(aes(group = Group), method = "wilcox.test",label = "p.signif") +
  theme(text = element_text(),
        axis.text.x = element_text(angle = 40,hjust = 1,vjust = 0.98),
        legend.text = element_text())
boxpplot
ggsave(plot = boxpplot,filename = "fig/step_3_miRNA_CIBERSORT.pdf",width = 8,height = 5)


# 保存数据
write.table(ssGSEAres,file = "output/step_3_ssGSEA_result.txt",
            row.names = T,col.names = T,quote = F,sep = "\t")


