##
##  标题：读取TCGA表达数据
##  作者：YHJ
##  时间：2021.11.3
##  备注：UCSC来源
## 

suppressMessages(library(stringr))
suppressMessages(library(dplyr))

Readdat <- function(inpath,idpath){
  print('Start reading and merging')
  ExpDat <- read.table(inpath,header = T,sep = '\t')
  colnames(ExpDat) <- gsub(colnames(ExpDat), pattern = '[.]', replacement = '_')
  names(ExpDat)[1] <- 'Ensembl_ID'
  symbol_id <- read.table(idpath,header = T,sep = '\t')
  names(symbol_id)[1] <- 'Ensembl_ID'
  ExpDat <- inner_join(symbol_id,ExpDat,by = 'Ensembl_ID')
  ExpDat <- aggregate(x = ExpDat[,7:ncol(ExpDat)],
                         by = list(ExpDat$gene),
                         FUN = median)
  rownames(ExpDat) <- ExpDat[,1]
  ExpDat <- ExpDat[,-1]
  ExpDat <- as.data.frame(ExpDat)
  print('successful')
  
  return(ExpDat)
}


ReadClinc <- function(inpath){
  print('Import clinical data')
  InClinc <- read.table(inpath,header = T,sep = '\t',quote = '')
  names(InClinc)[1] <- 'submitter'
  InClinc$submitter <- gsub(InClinc$submitter, pattern = '[-]', replacement = '_')
  InClinc$type <- ifelse(as.numeric(substr(InClinc$submitter,14,15)) > 9,'Normal','Tumor')
  rownames(InClinc) <- InClinc$submitter
  print('successful')
  
  return(InClinc)
}

# InClinc$OS.time <- signif(InClinc$days_to_death.demographic/30,3)
# InClinc$OS = ifelse(InClinc$vital_status.demographic == 'Alive',0,1)
# InClinc$Status <- ifelse(InClinc$OS == 1,'Dead','Alive')
# InClinc$Age <- InClinc$age_at_initial_pathologic_diagnosis
# InClinc$Gender <- InClinc$gender.demographic
# InClinc$Stage <- str_extract(InClinc$tumor_stage.diagnoses,pattern = "[iv]+")
# InClinc$lymphatic_invasion2 <- InClinc$lymphatic_invasion
# InClinc$nodal_deposits <- InClinc$non_nodal_tumor_deposits
# InClinc$diagnosis <- InClinc$number_of_first_degree_relatives_with_cancer_diagnosis
# InClinc$lymphnodes_he <- InClinc$number_of_lymphnodes_positive_by_he
# InClinc$lymphnod_ihc <- InClinc$number_of_lymphnodes_positive_by_ihc
