##
##  标题：GSEA分析
##  作者：YHJ
##  时间：2021.11.7
##  备注：只需要差异数据，可以选择本地或在线的
## 

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(enrichplot))


GSEaf <- function(indat,ingmt,inp,insep){
  GSEAdat <- subset(indat,P.Value < 0.05 & abs(logFC) > 1)
  kegmt <- read.gmt(ingmt)
  
  DEG_FC <- GSEAdat[,c('SYMBOL','logFC')]
  DEG_ENTIRE=bitr(DEG_FC$SYMBOL,fromType="SYMBOL",
                  toType="ENTREZID",
                  OrgDb="org.Hs.eg.db")
  DEG_ENTIRE <- DEG_ENTIRE[!duplicated(DEG_ENTIRE$SYMBOL),]
  GSEA_dat <- merge(DEG_ENTIRE,DEG_FC,by="SYMBOL")
  geneList <- GSEA_dat$logFC
  names(geneList) <- GSEA_dat$ENTREZID
  geneList <- sort(geneList,decreasing = T)
  # ## 在线绘制GSEA通路
  # GSEAA <- gseKEGG(geneList, organism = "hsa",pvalueCutoff = inp)
  ## 使用本地GMT文件
  GSEAA <- GSEA(geneList,pvalueCutoff = inp,TERM2GENE = kegmt)
  sortGSEAA <- GSEAA[order(GSEAA$enrichmentScore, decreasing = T),]
  write.csv(sortGSEAA,file = paste0('output/step_',insep,'_GSEA_Pthway.csv'))
  
  sortGSEAA <- sortGSEAA[sortGSEAA$pvalue < 0.05,]
  GSEA_PAY <- rownames(sortGSEAA)
  # GSEAplot0 <- gseaplot2(GSEAA, GSEA_PAY[1],
  #                        title = GSEA_PAY[1],
  #                        base_size = 8, ES_geom = "line",
  #                       subplots=c(1,2,3), pvalue_table = T,
  #                       color = "green")
  # GSEAplot0
  # pdf(file = paste0('fig/step_',insep,'_GSEA.pdf'),width = 6.5,height = 4)
  # print(GSEAplot0)
  # dev.off()
  
  for (i in 1:length(GSEA_PAY)) {
    gseap <- gseaplot2(GSEAA, GSEA_PAY[i],
                       title = GSEA_PAY[i],
                       base_size = 12, ES_geom = "line",
                       subplots=c(1,2,3), pvalue_table = T,
                       color = "green")
    
    pdf(file = paste0('fig/step_',insep,'_GSEA_', GSEA_PAY[i],'.pdf'),width = 8,height = 6)
    print(gseap)
    dev.off()
  }

}

