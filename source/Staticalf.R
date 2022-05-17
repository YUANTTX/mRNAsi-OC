##
##  标题：常用的建模分析
##  作者：YHJ
##  时间：2021.11.3
##  备注：uncox lasso mucox
##

suppressMessages(library(dplyr))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(glmnet))

uncoxf <- function(indat,ingene){
  un_mdat = indat
  covariates <- gsub(ingene, pattern = '[-]', replacement = '_')
  colnames(un_mdat) <- gsub(colnames(un_mdat), pattern = '[-]', replacement = '_')
  
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = un_mdat)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value   <- signif(x$wald["pvalue"], digits=3)
                           wald.test <- signif(x$wald["test"], digits=3)
                           beta      <- signif(x$coef[1], digits=3)
                           HR        <- sprintf("%0.3f", x$coef[2])
                           HR95L     <- sprintf("%0.3f", x$conf.int[,"lower .95"])
                           HR95H     <- sprintf("%0.3f", x$conf.int[,"upper .95"])
                           HRCI      <- paste0(HR, "(", HR95L, "-", HR95H, ")")
                           res       <- c(beta, wald.test,HR,HR95L,HR95H, p.value,HRCI)
                           names(res)<- c("Beta","wald.test",'HR','HR95L','HR95H', "pvalue","HR_95_CI")
                           return(res)
                         })
  uncox <- as.data.frame(univ_results)
  uncox <- t(as.data.frame(univ_results, check.names = FALSE))
  uncox <- as.data.frame(uncox)
  rownames(uncox) <- gsub(rownames(uncox), pattern = '_', replacement = '-')
  uncox$pvalue <- as.numeric(uncox$pvalue)
  
  return(uncox)
}

lassof <- function(indat,ingene,insep){
  inGene = ingene
  lasso_dat = indat
  colnames(lasso_dat) <- gsub(colnames(lasso_dat), pattern = '[-]', replacement = '_')
  x <- as.matrix(lasso_dat[,gsub(inGene, pattern = '-', replacement = '_')])
  y <- lasso_dat[,c('OS.time', 'OS')]
  y$OS.time[y$OS.time<0.1] = 0.1
  y$OS.time <- as.double(y$OS.time)
  y$OS <- as.double(y$OS)
  y <- y[!is.na(y$OS.time),]
  x <- x[rownames(y),]
  
  set.seed(1)
  fit <- glmnet(x,Surv(y$OS.time, y$OS),family = 'cox',alpha = 1, nlambda = 100)
  pdf(file = paste0('fig/step_',insep,'_lasso1.pdf'),width = 4,height = 4)
  plot(fit,label = T)
  dev.off()
  
  set.seed(1)
  y <- as.matrix(survival::Surv(y$OS.time, y$OS))
  # lasso_fit <- cv.glmnet(x, y, family='cox', type.measure = 'deviance',nfolds = nrow(x))
  lasso_fit <- cv.glmnet(x, y, family='cox', type.measure = 'deviance',nfolds = 100)
  pdf(file = paste0('fig/step_',insep,'_lasso2.pdf'),width = 4,height = 4)
  plot(lasso_fit)
  dev.off()
  coefficient <- coef(lasso_fit, s=lasso_fit$lambda.min)
  Active.Index <- which(as.numeric(coefficient) != 0)
  lassodat = data.frame(coef = as.numeric(coefficient)[Active.Index],
                        gene = rownames(coefficient)[Active.Index])
  lassodat$gene <- gsub(lassodat$gene, pattern = '_', replacement = '-')
  rownames(lassodat) = lassodat$gene
  return(lassodat)
}

mucoxf <- function(indat,ingene){
  mu_mdat =   indat
  signagene = gsub(ingene, pattern = '[-]', replacement = '_')
  colnames(mu_mdat) <- gsub(colnames(mu_mdat), pattern = '[-]', replacement = '_')
  multivariate <- as.formula(paste0('Surv(OS.time,OS)~', 
                                    paste(signagene, sep = '', collapse = '+')))
  res.cox <- coxph(multivariate, data = mu_mdat)
  multiCoxSum <- summary(res.cox)
  out_multi <- data.frame()
  out_multi <- cbind(
    coef     = multiCoxSum$coefficients[,"coef"],
    HR       = multiCoxSum$conf.int[,"exp(coef)"],
    HR95L   = multiCoxSum$conf.int[,"lower .95"],
    HR95H   = multiCoxSum$conf.int[,"upper .95"],
    pvalue   = multiCoxSum$coefficients[,"Pr(>|z|)"])
  out_multi <- as.data.frame(out_multi)
  out_multi$pvalue <- signif(out_multi$pvalue, digits=3)
  out_multi$HR = sprintf("%0.3f", out_multi$HR)
  out_multi$HR95L = sprintf("%0.3f", out_multi$HR95L)
  out_multi$HR95H = sprintf("%0.3f", out_multi$HR95H)
  out_multi$HR_95_CI = paste0(out_multi$HR, "(", out_multi$HR95L, "-", out_multi$HR95H, ")")

  mucox <- as.data.frame(out_multi)
  mucox$pvalue <- as.numeric(mucox$pvalue)
  rownames(mucox) <- gsub(rownames(mucox), pattern = '_', replacement = '-')
  return(mucox)
}