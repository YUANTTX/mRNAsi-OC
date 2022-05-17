##
##  标题：使用本地注释转换ID
##  作者：yhj
##  时间：2021.9.2
##

suppressMessages(library(dplyr))
suppressMessages(library(stringr))

Aff2Symbol <- function(indat,iddat,fromID,toID){
  # AffID <- read.table('GPL570-55999.txt',header = T,fill = T,quote = '',sep = '\t',comment.char = '#')
  ExpDat <- as.data.frame(indat)
  ExpDat$name <- rownames(ExpDat)
  names(ExpDat)[names(ExpDat) == 'name'] <- fromID
  AffID <- iddat[,c(fromID,toID)]
  AffID[,1] <- str_split(AffID[,1],'[ /// ]',simplify = T)[,1]
  AffID[,2] <- str_split(AffID[,2],'[ /// ]',simplify = T)[,1]
  
  ExpDat <- inner_join(AffID,ExpDat,by = fromID)
  ExpDat <- ExpDat[!duplicated(ExpDat[,toID]),]
  rownames(ExpDat) <- ExpDat[,toID]
  ExpDat <- ExpDat[,-c(1:2)]

  return(ExpDat)
}

# Ensembl2symbol <- function(inpath,idpath){
#   suppressMessages(library(stringr))
#   
#   print('Converting to gene name')
#   Expdat <- read.table(inpath,header = T,sep = '\t')
#   Symbol <- read.table(idpath,header = T,sep = '\t')
#   names(Symbol)[1] <- 'Ensembl_ID'
#   Expdat <- inner_join(Symbol,Expdat,by = 'Ensembl_ID')
#   Expdat <- Expdat[!duplicated(Expdat$gene),]
#   rownames(Expdat) <- Expdat$gene
#   Expdat <- Expdat[,-c(1,3:6)]
#   names(Expdat)[1] <- 'gene_id'
# 
#   print('successful')
#   return(Expdat)
# }

## 附录（使用方法）
# Gencode <- read.table('/home/yuan/DataBase/UCSC/GencodeID/gencode.v22.annotation.gene.probeMap',
#                       header = T,sep = '\t')
# TCGX_EXP <- Aff2Symbol(indat = TCGX_EXP,iddat = Gencode,fromID = 'id',toID = 'gene')