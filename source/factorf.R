##
##  标题：风险因子关联图
##  作者：YHJ
##  时间：2021.11.3
##

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(ggplotify))
suppressMessages(library(cowplot))
suppressMessages(library(Hmisc))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(gridExtra))
suppressMessages(library(circlize))


factorf <- function(indat,inGene,insep){
  
  indat <- indat[order(indat$RiskScore),]
  s = as.formula(paste('Surv(OS.time, OS)~', noquote(paste(inGene,collapse = ' + '))))
  model <- coxph(s, data = indat )
  summary(model,data=indat)
  RiskScore <- indat$RiskScore
  names(RiskScore) = rownames(indat)
  fp <- RiskScore
  phe <- indat
  fp_dat = data.frame(patientid = 1:length(fp),
                      fp        = as.numeric(sort(fp)))
  fp_dat$riskgroup= ifelse(fp_dat$fp >= median(fp_dat$fp),'high','low')
  
  sur_dat=data.frame(patientid = 1:length(fp),
                     time      = phe[names(sort(fp)),'OS.time'],
                     status    = phe[names(sort(fp )),'OS'])
  sur_dat$status = ifelse(sur_dat$status==0,'Alive','Dead')
  sur_dat$status = factor(sur_dat$status,levels = c("Dead","Alive"))
  sur_dat$time <- signif(sur_dat$time/12,3)
  exp_dat <- indat[names(sort(fp)),
                     c(inGene)]
  p1 <- ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup),size = 1)+
    scale_colour_manual(values = c("#F71E35","#00a03e"))+
    theme_test()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
    geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),
               colour="black", linetype="dotted",size=0.8) +
    theme(axis.title.x=element_blank(),axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
          legend.text=element_text(size=14),legend.title =element_text(size=14))
  ggsave(plot = p1,filename = paste0('fig/step_',insep,'_fator_p1.pdf'),width = 6,height = 3)
  p2 <- ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=status),size = 1.5)+
    theme_test()+
    scale_colour_manual(values = c("#F71E35","#63b931"))+
    labs(x="Patient ID(increasing risk score)",y="Survival time(years)")+
    geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", 
               linetype="dotted",size=0.8) +
    theme(axis.title.x=element_blank(),axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
          legend.text=element_text(size=14),legend.title =element_text(size=14))
  ggsave(plot = p2,filename = paste0('fig/step_',insep,'_fator_p2.pdf'),width = 6,height = 3)
  # mycolors <- colorRampPalette(c("#304ffe", "white", "#FF6600"), bias = 1.2)(100)
  # mycolors <- colorRamp2(c(-1.5, 0, 1.5), c("#99CC33", "white", "#FF9900"))
  mycolors <- colorRamp2(c(-2, 0, 2),c("#99CC33", "white", "red"))
  tmp=t(scale(exp_dat))
  tmp[tmp > 2] = 2
  tmp[tmp < -2] = -2
  p3 = pheatmap(tmp,name = ' ',
                cluster_cols = F,cluster_rows = T,
                show_rownames = T,show_colnames = F,
                treeheight_row = 12,
                col= mycolors, border_color = F)
  heig = signif(length(rownames(tmp))/3,3)
  pdf(file = paste0('fig/step_',insep,'_fator_p3.pdf'),width = 6,height = heig)
  print(p3)
  dev.off()
  
  pdf(file = paste0('fig/step_',insep,'_fator_all.pdf'),width = 8,height = 10)
  plots = list(p1,p2,as.ggplot(as.grob(p3)))
  lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7)))
  grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
  # plots = list(p1,as.ggplot(as.grob(p3)))
  # lay1 = rbind(c(rep(1,7)),c(rep(2,7)))
  # grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10))
  dev.off()
}


