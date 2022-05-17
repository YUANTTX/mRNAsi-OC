##
##  标题：基因差异分析
##  作者：yhj
##  时间：2021.9.2
##


limmaf <- function(indat,ingroup,inp,inlogfc,batoh){
  
  suppressMessages(library(limma))
  suppressMessages(library(dplyr))
  
  print('limma is being used for differential analysis')
  design <- model.matrix(~0+factor(ingroup))
  colnames(design) = levels(factor(ingroup))
  rownames(design) = colnames(indat)
  con.matrix <- makeContrasts(batoh,levels = design)
  fit <- lmFit(indat,design)
  fit2 <- contrasts.fit(fit,con.matrix)
  fit2 <- eBayes(fit2)
  temput <- topTable(fit2,coef = 1,n = Inf)
  temput[sapply(temput,is.infinite)] <- NA
  temput <- na.omit(temput)
  temput$State <- ifelse(temput$adj.P.Val > inp,'not',
                         ifelse(temput$logFC > inlogfc,'up',
                                ifelse(temput$logFC < -inlogfc,'down','not')))
  temput$SYMBOL <- rownames(temput)
  print('successful')
  return(temput)
}




DEseq2f <- function(indat,ingroup,inp,inlgfc){
  
  suppressMessages(library(dplyr))
  suppressMessages(library(DESeq2))
  
  countData <- indat
  condition <- factor(ingroup)
  ids <- rownames(countData)
  countData <- apply(as.matrix(countData),2,function(x) as.integer(x))
  rownames(countData) <- ids
  countData <- as.data.frame(countData)
  
  print('DESeq2 is being used for differential analysis')
  colData <- data.frame(row.names=colnames(countData), condition)
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
  dds <- DESeq(dds)
  res <- results(dds)
  resinfo <- as.data.frame(res)
  resinfo$SYMBOL <- rownames(resinfo)
  resinfo$State <- ifelse(resinfo$pvalue > inp,'not',
                          ifelse(resinfo$log2FoldChange >= inlgfc,'up',
                                 ifelse(resinfo$log2FoldChange <= -inlgfc,'down','not')))
  names(resinfo)[names(resinfo) == 'log2FoldChange'] <- 'logFC'
  names(resinfo)[names(resinfo) == 'pvalue'] <- 'P.Value'
  print('successful')
  return(resinfo)
}



Edgrf <- function(indat,ingroup,inp,inlgfc){
  
  suppressMessages(library(dplyr))
  suppressMessages(library(edgeR))
  
  DEG <- DGEList(counts=indat,group=factor(ingroup))
  keep <- rowSums(cpm(DEG)>1) >= 2
  DEG <- DEG[keep, , keep.lib.sizes=FALSE]
  DEG$samples$lib.size <- colSums(DEG$counts)
  DEG <- calcNormFactors(DEG)
  design <- model.matrix(~0+factor(ingroup))
  rownames(design) <- colnames(DEG)
  colnames(design) <- levels(factor(ingroup))
  DEG <- estimateGLMCommonDisp(DEG,design)
  DEG <- estimateGLMTrendedDisp(DEG, design)
  DEG <- estimateGLMTagwiseDisp(DEG, design)
  fit <- glmFit(DEG, design)
  lrt <- glmLRT(fit, contrast=c(1,-1)) 
  EDGR_INFO <- topTags(lrt, n=nrow(DEG))
  EDGR_INFO <- as.data.frame(EDGR_INFO)
  EDGR_INFO$Gene <- rownames(EDGR_INFO)
  EDGR_INFO$State <- ifelse(EDGR_INFO$PValue > inp,'not',
                            ifelse(EDGR_INFO$logFC >= inlgfc,'up',
                                   ifelse(EDGR_INFO$logFC <= -inlgfc,'down','not')))
  return(EDGR_INFO)
}
