##
##    step7 新亚型与免疫关系
## 

## 空间设置
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(library(GSEABase))
suppressMessages(library(ggstatsplot))
suppressMessages(library(rlang))
source('source/immunef.R')
source('source/DEGplotf.R')
source('ImmDAT/CIBERSORT/Cibersort.R')


## 加载数据
load("rdata/step_1_TCGA_EXP.Rdata")
load('rdata/step_5_TCGA_clinc.Rdata')
load("ImmDAT/imm1316/ImmGene_Set.Rdata")
geneSets <- getGmt("ImmDAT/imm1316/Hhimmune.symbol.gmt")
checkGene <- read.csv('ImmDAT/ImmCheckGene.csv')

## 整理数据
TCGA_clinc <- TCGA_clinc[order(TCGA_clinc$StemType),]
TCGA_EXP <- TCGA_EXP[,TCGA_clinc$submitter]

##  亚型与ESTIMATE得分
names(TCGA_clinc)
color2 = c('#e19d29','#4298b5')
comparison = list(c("Steamness Subtype I","Steamness Subtype II"))

vars = c('ImmuneScore','StromalScore','ESTIMATEScore')
for (i in vars) {
  boxp9 = ggboxplot(TCGA_clinc, size = 1, bxp.errorbar = F,
                    x="StemType", y= i, color = "StemType",
                    palette = color2,
                    order = c("Steamness Subtype I","Steamness Subtype II")) + 
    stat_compare_means(comparisons = comparison, method = "t.test") +
    theme(text = element_text(size = 10),legend.text = element_text(size = 14),
          legend.position = 'none') + 
    scale_x_discrete('StemType',labels = c("Steamness Subtype I" = "Subtype I","Steamness Subtype II" = "Subtype II"))
  print(boxp9)
  pdf(file=paste0("fig/step_7_Stemtype_",i,".pdf"),width = 3,height = 4)
  print(boxp9)
  dev.off()
}


## 亚型与22种免疫细胞相关性
LM22.file <- "ImmDAT/CIBERSORT/LM22.txt"
TCGA_exp.file = 'tmp/step_1_TCGA_EXP.txt'

TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 1000, QN = F)
TCGA_TME.results <- TCGA_TME.results[TCGA_clinc$submitter,]

group_list <- TCGA_clinc$StemType %>% factor(.,levels = c("Steamness Subtype I","Steamness Subtype II"))

TME_data <- as.data.frame(TCGA_TME.results[,1:22])
TME_data$group <- group_list
TME_data$sample <- row.names(TME_data)
TME_New = melt(TME_data)
colnames(TME_New)=c("Group","Sample","Celltype","Immune cell composition")

boxpplot <- ggboxplot(TME_New, x = "Celltype", y = "Immune cell composition",
                      fill = 'Group', palette = c('#e19d29','#4298b5'),
                      xlab = F) + 
  stat_compare_means(aes(group = Group), method = "wilcox.test",label = "p.signif") +
  theme(text = element_text(),
        axis.text.x = element_text(angle = 40,hjust = 1,vjust = 0.98),
        legend.text = element_text())
boxpplot
ggsave(plot = boxpplot,filename = "fig/step_7_SubType_CIBERSORT.pdf",width = 8,height = 6)


## 免疫检查点表达 
TCGA_ck <- as.data.frame(t(TCGA_EXP[checkGene$symbol,]))
TCGA_ck$submitter <- rownames(TCGA_ck)
TCGA_clinc2 <- inner_join(TCGA_clinc,TCGA_ck,by = "submitter")
names(TCGA_clinc2)
TCGA_clinc2$StemType

color2 = c('#e19d29','#4298b5')
comparison = list(c("Steamness Subtype I","Steamness Subtype II"))
for (i in c("PDCD1","CD274","PDCD1LG2","CTLA4","CD80","CD86")){
  boxp2 = ggboxplot(TCGA_clinc2, size = 1, bxp.errorbar = F,
                    x="StemType", y = i,
                    title = i,ylab = 'Expression (log2(FPKM+1))',xlab = F,
                    fill = 'StemType', palette = c('#0099CC','#996699'),
                    order = c("Steamness Subtype I","Steamness Subtype II")) + 
    stat_compare_means(comparisons = comparison, method = "wilcox.test") +
    theme(text = element_text(size = 10),legend.text = element_text(size = 14),
          legend.position = 'none') + 
    scale_x_discrete('StemType',labels = c("Steamness Subtype I" = "Subtype I","Steamness Subtype II" = "Subtype II"))
  print(boxp2)
  pdf(file = paste0('fig/step_7_StemType_',i,'.pdf'),width = 2.6,height = 3.6)
  print(boxp2)
  dev.off()
}
