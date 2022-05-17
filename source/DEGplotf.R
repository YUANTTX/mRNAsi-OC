##
##  标题：差异分析结果可视化
##  作者：YHJ
##  时间：2021.11.3
##

suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggplot2))


DEGPlotf <-  function(indat,ininfo,intype,ingroup,insep){
  DEGdat <- indat
  DEGinfo <- ininfo
  grouplist <- ingroup
  
  print('Drawing volcano map and heat map')
  DEGinfo$P.value2 = -log10(DEGinfo$P.Value)
  plot1 <- ggplot(DEGinfo, aes(x = logFC, y = P.value2, colour=State)) +
    geom_point(alpha=0.7, size=1) + 
    # scale_color_manual(values=c("#546de5", "black","#ff4757"))+
    scale_color_manual(values=c("#00a03e", "black","#ff4757"))+
    geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
    labs(x="logFC",y="-log10(Pvalue)")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank(),
          legend.key = element_rect(colour ="black"))
  plot1
  ggsave(plot = plot1,filename = paste0('fig/step_',insep,'_DEG_vol.pdf'),width = 6,height = 6)
  
  FC = DEGinfo$logFC
  names(FC) = rownames(DEGinfo)
  FCgene = c(names(head(sort(FC),50)),names(tail(sort(FC),50)))
  ht_dat = indat[FCgene,]
  ann_col = data.frame(Groups=grouplist)
  rownames(ann_col) = colnames(ht_dat)
  Groupscol <- c('#248888','#E7475E')
  names(Groupscol) <- intype
  ann_colors = list(Groups = Groupscol)
  main_col = colorRamp2(c(-3, 0, 3), c("#00e701", "black", "#eb0c08"))
  plot2 <- pheatmap(as.matrix(ht_dat),name = ' ',
                    show_rownames = F,show_colnames = F,
                    cluster_rows= T,cluster_cols = F,
                    treeheight_col = 10,
                    # main="Heatmap of top 200 lncRNA",
                    annotation_row = NA,
                    annotation_col = ann_col,
                    annotation_colors =ann_colors,
                    border_color = F,
                    scale="row",color = main_col)
  plot2
  pdf(file = paste0('fig/step_',insep,'_DEG_heatmap.pdf'),width = 8,height =7)
  print(plot2)
  dev.off()
  print('successful')
}