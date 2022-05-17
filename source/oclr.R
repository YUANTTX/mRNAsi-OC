##
##  标题：将干性指数计算封装到函数中
##  作者：YHJ
##  时间：2021.11.3
##

suppressMessages(library(gelnet))
suppressMessages(library(dplyr))
suppressMessages(library(biomaRt))
suppressMessages(library(synapser))



genes2hugo <- function(v, srcType = "ensembl_gene_id"){
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl <-biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             host="www.ensembl.org",
                             dataset="hsapiens_gene_ensembl")
  ID <-biomaRt::getBM(attributes=c(srcType, "hgnc_symbol"),
                      filters=srcType,
                      values=v, mart=ensembl)
  ##Make sure there was at least one mapping
  if(nrow(ID) < 1) top("No IDs mapped successfully")
  ## Drop empty duds
  j <- which(ID[,2] == '')
  if(length(j) >0 ) ID <- ID[-j,]
  stopifnot(all(ID[,1] %in% v))
  ID
}


f <- function(v){
  unlist(lapply(strsplit(v,'\\|'),'[[',1))
}

Downpcbc <- function(oupath){
  
  print('Download and calculate dryness index model')
  synLogin("yuanttx","ubuntu22")
  synRNA <- synGet("syn2701943",downloadLocation = paste0(oupath,'/pcbc'))
  X <- read.delim(synRNA$path) %>%
    tibble::column_to_rownames('tracking_id') %>%
    as.matrix()
  synMate <- synTableQuery('SELECT UID, Diffname_short FROM syn3156503')
  synMate <- as.data.frame(synMate)
  Y <- synMate %>%
    mutate(UID = gsub('-','.',UID)) %>%
    tibble::column_to_rownames('UID')
  y = Y[colnames(X),]
  y['SC11.014beb.133.5.6.11'] <- 'EB'
  y['SC12.039ECTO.420.436.92.16'] <- 'ECTO'
  v <- strsplit(rownames(X),'\\.') %>% lapply('[[',1) %>% unlist()
  rownames(X) <- v
  V <- genes2hugo(rownames(X))
  X <- X[V[,1],]
  rownames(X) <- V[,2]
  fnGenes = NULL
  if(is.null(fnGenes) == FALSE){
    vGenes <- read.delim(fnGenes,header = F) %>% as.matrix() %>% drop()
    VE <- genes2hugo(vGenes,'entrezgene')
    X <- X[intersect(rownames(X),VE[,2]),]
  }
  m <- apply(X, 1, mean)
  X <- X -m
  j <- y[y$Diffname_short == 'SC',]
  X.tr <- X[,colnames(X) %in% rownames(j)]
  X.bk <- X[,!colnames(X) %in% rownames(j)]
  mm <- gelnet(t(X.tr),NULL,0,1)
  w = mm$w
  write.table(mm$w,file = paste0(oupath,'pcbc-stemsig.tsv'),sep = '\t',quote = F,col.names = F)
  print('successful')
}


oclrpd <- function(sigpath,datpath){
  w <- read.delim(sigpath,header = F,row.names = 1) %>%
    as.matrix() %>%
    drop()
  X <- read.delim(datpath,as.is = T,check.names = F) %>%
    filter(!grepl('\\?',gene_id)) %>%
    mutate(gene_id = f(gene_id)) %>%
    filter(gene_id %in% names(w))
  j <- grep('SCL35E2',X[,1])
  if(length(j) > 1)
    X <- X[-j[-1],]
  rownames(X) <- NULL
  X <- X %>% tibble::column_to_rownames('gene_id') %>% as.matrix()
  stopifnot(all(rownames(X) %in% names(w)))
  w <- w[rownames(X)]
  s <- apply(X,2,function(z) {cor(z,w,method = 'sp',use = 'complete.obs')})
  s <- s - min(s)
  s <- s / max(s)
  
  redat <- as.data.frame(cbind(s))
  names(redat)[1] <- 'StenScore'
  redat$StenGroup <- ifelse(redat$StenScore > median(redat$StenScore),'high','low')
  redat$mRNAsi <- redat$StenScore
  redat$mRNAsiGroup <- redat$StenGroup
  redat$submitter <- rownames(redat)
  return(redat)
}



