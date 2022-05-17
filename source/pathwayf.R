##
##  标题：通路富集分析：GO & KEGG
##  作者：YHJ
##  时间：2021.11.3
##  备注：注意需不需要区分高低表达分组；(需不需要差异信息)
##

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

Pathwayf <- function(ingene,insep){
  
  payinfo <- data.frame('SYMBOL' = ingene)
  print('Converting ID ...')
  nrDEG <- bitr(unique(payinfo$SYMBOL),  OrgDb = org.Hs.eg.db,
                fromType = "SYMBOL",toType = c("ENTREZID"))
  nrDEG <- merge(payinfo,nrDEG,by.x = 'SYMBOL',by.y = 'SYMBOL')
  # gene_up   <- nrDEG[nrDEG$State == 'up','ENTREZID']
  # gene_down <- nrDEG[nrDEG$State == 'down','ENTREZID']
  gene_all  <- as.character(nrDEG[,'ENTREZID'])
  
  print('Analyzing pathway enrichment')
  tryout <- try(   go.all <- enrichGO(gene          = gene_all, 
                                      OrgDb         = org.Hs.eg.db,
                                      ont           = 'all', 
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.5, 
                                      qvalueCutoff  = 0.5, 
                                      readable      = TRUE),
                silent = T)
  if("try-error" %in% class(tryout)){
    print('There is a problem, please try again.')
  }else{
    goall <- barplot(go.all,split="ONTOLOGY",title = 'Enrichment GO') + 
      facet_grid(ONTOLOGY~., scale="free") +
      scale_y_discrete(labels=function(x) str_wrap(x, width=60))
    pdf(file = paste0('fig/step_',insep,'_GO.pdf'),width = 8,height = 8)
    print(goall)
    dev.off()
  }
  
  tryout <- try(  kk.all <- enrichKEGG(gene         = gene_all,
                                       organism     = 'hsa',
                                       pvalueCutoff = 0.5,
                                       qvalueCutoff = 0.5),
                silent = T)
  if("try-error" %in% class(tryout)){
    print('There is a problem, please try again.')
  }else{
    kkall <- barplot(kk.all,showCategory=30, title ="Enrichment KEGG") +
      scale_y_discrete(labels=function(x) str_wrap(x, width=60))
    pdf(file = paste0('fig/step_',insep,'_KEGG.pdf'),width = 8,height = 8)
    print(kkall)
    dev.off()
  }
  
  tryout <- try(  save(go.all,kk.all, file = paste0('rdata/step_',insep,'_pathway.Rdata')),
                  silent = T)
  if("try-error" %in% class(tryout)){
    print('There is a problem, please try again.')
  }else{
    print('successful')
  }
}


DavidTabf <- function(indat){
  DaviDTab <- indat
  DaviDTab <- separate(DaviDTab, Term, sep = "~",into = c("ID", "Term"))
  DaviDTab <- DaviDTab[,c('Term','Genes')]
  
  DaviDTab$Genes <- gsub(DaviDTab$Genes, pattern = '"', replacement = '')
  GenesTab = str_split(DaviDTab$Genes ,', ',simplify = T)
  OutTab <- data.frame()
  for (i in 1:nrow(DaviDTab)) {
    TableD <- data.frame(Term = DaviDTab[i,1], Genes = GenesTab[i,])
    OutTab <- rbind(OutTab,TableD)
  }
  OutTab <- OutTab[OutTab$Genes != '',]
  
  return(OutTab)
}


DavidPlotf <- function(indat){
  keggSig = keggSig[keggSig$PValue < 0.05,]
  
  ## 2.把Term拆分为通路与描述c("ID", "Term")
  library(tidyr)
  keggSig = separate(keggSig, Term, sep = ":",
                     into = c("ID", "Term"))
  
  library(ggplot2)
  ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
    geom_point(aes(size=Count,color=-1*log10(PValue)))+
    scale_colour_gradient(low="green",high="red")+
    labs(
      color=expression(-log[10](P.value)),
      size="Gene number",
      x="Fold enrichment"
      # y="Pathway name",
      # title="Pathway enrichment")
    )+
    theme_bw()+
    theme(
      axis.text.y = element_text(size = rel(1.3)),
      axis.title.x = element_text(size=rel(1.3)),
      axis.title.y = element_blank()
    )
  ggsave('fig/step_2_kegg.pdf',width = 7,height = 8)
  
  ## 上调基因
  up=read.table('input/david.txt',sep = '\t',header = TRUE,quote = '')
  up_rt=up[up$PValue< 0.05,]
  library(tidyr)
  up_rt = separate(up_rt, Term, sep = "~",
                   into = c("ID", "Term"))
  ## 上调基因的bp
  bp_df = up_rt[up_rt$Category == 'GOTERM_BP_DIRECT',]
  bp_df = bp_df[order(bp_df$Count,decreasing = T),]
  bp = bp_df[1:10,]
  ## 上调基因的cc
  cc_df = up_rt[up_rt$Category == 'GOTERM_CC_DIRECT',]
  cc_df = cc_df[order(cc_df$Count,decreasing = T),]
  cc = cc_df[1:10,]
  ## 上调基因的mf
  mf_df = up_rt[up_rt$Category == 'GOTERM_MF_DIRECT',]
  mf_df = mf_df[order(mf_df$Count,decreasing = T),]
  mf = mf_df[1:10,]
  ## 上调基因的all
  allGo = rbind(bp,cc,mf)
  library(stringr)
  table(allGo$Category)
  allGo$Category = substr(allGo$Category,8,9)
  
  ## 绘制上调基因三种类别
  library(ggpubr)
  colnames(allGo)
  p = ggbarplot(data = allGo,x = "ID",y = 'Count',
                fill = "Category",
                palette = c("cadetblue3","mediumslateblue","mediumorchid3"),
                sort.by.groups = T,xlab = '',ylab = "Target genes",x.text.angle = 60)
  p
  ggsave(plot = p,'barplot.pdf',width = 10,height = 10)
}