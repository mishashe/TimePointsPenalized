#' TimePointsPenalized: Lasso with penalized differences between adjacent time points 
#'
#' 
#' @section Mypackage functions:
#' A package to fit lasso with penalized differences between adjacent time points coefficients.
#'
#' @docType package
#' @name TimePointsPenalized
#' @useDynLib TimePointsPenalized, .registration=TRUE
NULL



#' Fit lasso with penalized differences between adjacent time points coefficients 
#'
#' @param y0 case/control vector (no time iformation - "naive" approach)
#' @param x0 gene expression matrix (rows-samples by columns-genes)
#' @param FollowUp follow-up times (recurrence time for recurrences and follow-up for patients with no recurrences)
#' @param lam1V array of lasso penalty prefactor
#' @param gamma prefactor of the second penalty term - differences between adjacent time points coefficients
#' @param tV array of time points
#' @param Clinical0 dataframe with clinical information (same order as rows of x0)
#' @export
fitTimePointsPenalized <- function(y0, x0, FollowUp, lam1V, gamma, tV, 
                                   Clinical0=data.frame(case_control0=y0), startWithGlmnet=FALSE){  
  if (startWithGlmnet){
    fits0 <- fitTimePointsNonPenalized(y0, x0, FollowUp, lam1V, gamma, tV, Clinical0=data.frame(case_control0=y0))
  }
  else{
    fits0 <-  NULL
  }
  Intercept <- rep(0,length(tV))
  beta <- rep(0,ncol(x0)*length(tV))
  y <- c()
  Clinical <- data.frame()
  samplesT <- 1:(nrow(x0)*length(tV))
  GenesT <- rep("",ncol(x0)*length(tV))
  Clinical0$sample <- rownames(x0)
  Clinical0$FollowUp <- FollowUp
  
  # make long Clinical and y vectors containing all time points
  for (it in 1:length(tV)){
    t <- tV[it]
    case_controlT <- ifelse(FollowUp>t,0,ifelse(y0==1,1,-1))
    if (sum(case_controlT==0)<2){
      print(paste0("Not enough controls for t=",t))
      print(table(case_controlT))
      return(NULL)
    }
    if (sum(case_controlT==1)<2){
      print(paste0("Not enough cases for t=",t))
      print(table(case_controlT))
      return(NULL)
    }
    y <- c(y,case_controlT)
    ClinicalT <- Clinical0
    ClinicalT$StatusT <- case_controlT
    ClinicalT$FollowUp <- FollowUp
    ClinicalT$time <- tV[it]
    ClinicalT$samples <- rownames(x0)
    ClinicalT$samplesT <- paste0(rownames(x0),"_t_",it)
    Clinical <- rbind(Clinical,ClinicalT)
    GenesT[(1+(it-1)*ncol(x0)):(it*ncol(x0))] <- paste0(colnames(x0),"_t_",it)
    samplesT[(1+(it-1)*nrow(x0)):(it*nrow(x0))] <- paste0(rownames(x0),"_t_",it)
  }
  names(beta) <- GenesT
  names(y) <- samplesT
  Ind <- which(y %in% c(0,1))
  Clinical <- Clinical[Ind,]
  y <- y[Ind]
  
  IndFor0 <- c()
  IndTFor0 <- c()
  w <- y*0
  for (it in 1:length(tV)){
    IndT <- which(Clinical$time==tV[it])
    w[IndT][y[IndT]==0] <- 1/sum(y[IndT]==0) #weights balance the total weight of cases and controls
    w[IndT][y[IndT]==1] <- 1/sum(y[IndT]==1)
    w[IndT] <- w[IndT]/sum(w[IndT])*length(IndT)
    IndFor0 <- c(IndFor0,which(rownames(x0) %in% Clinical$samples[IndT]))
    IndTFor0 <- c(IndTFor0,which(rownames(x0) %in% Clinical$samples[IndT])*0+it)
  }
  w <- w/sum(w)
  fits <- list(list())
  for (it in 1:length(tV)){
    fits[[it]] <- list()
    fits[[it]]$beta <- matrix(0,nrow=dim(x0)[2],ncol=0)
    fits[[it]]$Intercept <- c()
    fits[[it]]$lambda <- c()
    fits[[it]]$gamma <- c()
  }
  for (ilam1 in 1:length(lam1V)){
    lam1 <- lam1V[ilam1]/length(tV)
    lam2 <- gamma*lam1
    if (!is.null(fits0)){
      for (it in 1:length(tV)){
        beta[(1:(dim(x0)[2]))+(it-1)*dim(x0)[2]] <- fits0[[it]]$beta[,ilam1,drop=TRUE]
        Intercept[it] <- fits0[[it]]$a0[ilam1]
      }
    }
    fit <- Fit(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0, IndTFor0)
    beta <- fit$beta
    Intercept <- fit$Intercept
    for (it in 1:length(tV)){
      betaOut <- fit$beta[(1:(dim(x0)[2]))+(it-1)*dim(x0)[2],1,drop=TRUE]
      InterceptOut <- fit$Intercept[it]
      IndT <- which(Clinical$time==tV[it])
      fits[[it]]$beta <- cbind(fits[[it]]$beta,betaOut)
      fits[[it]]$Intercept <- c(fits[[it]]$Intercept,InterceptOut)
      fits[[it]]$lambda <- c(fits[[it]]$lambda,lam1)
      fits[[it]]$gamma <- c(fits[[it]]$gamma,gamma)
    }
  }
  return(fits)
}       

#' CV with penalized differences between adjacent time points coefficients 
#'
#' @param y0 case/control vector (no time iformation - "naive" approach)
#' @param x0 gene expression matrix (rows-samples by columns-genes)
#' @param FollowUp follow-up times (recurrence time for recurrences and follow-up for patients with no recurrences)
#' @param lam1V array of lasso penalty prefactor
#' @param gamma prefactor of the second penalty term - differences between adjacent time points coefficients
#' @param tV array of time points
#' @param Clinical0 dataframe with clinical information (same order as rows of x0)
#' @export
fitTimePointsPenalized.cv <- function(y0, x0, FollowUp, lam1V, gamma, tV, Clinical0=data.frame(case_control0=y0), 
                                      startWithGlmnet=FALSE,folds,whatToMaximize="auc"){
  dataCV <- foreach (fold = c(-1,unique(folds)), .inorder=FALSE) %dopar%
  {
    if (fold==-1){
      return(fitTimePointsPenalized(y0, x0, FollowUp, lam1V, gamma, tV, Clinical0=data.frame(case_control0=y0), startWithGlmnet))
    }
    print(fold)
    Ind <- which(fold!=folds)
    fits <- fitTimePointsPenalized(y0[Ind], x0[Ind,], FollowUp[Ind], lam1V, gamma, tV, Clinical0=data.frame(case_control0=y0[Ind]), startWithGlmnet)
    data <- data.frame()
    for (it in 1:length(tV))
    {
      status <- ifelse(FollowUp[-Ind] > tV[it],0,ifelse(y0[-Ind]==1,1,-1))
      dataT <- data.frame(sample=rownames(x0)[-Ind], status=status,type=y0[-Ind], timepoint=tV[it], FollowUp=FollowUp[-Ind])
      for (ilam1 in 1:length(lam1V))
      {
        beta <- fits[[it]]$beta[,ilam1,drop=FALSE]
        Intercept <- fits[[it]]$Intercept[ilam1]
        preds <- round(1/(1+exp(-x0[-Ind,] %*% beta - Intercept)),2)
        preds[preds<1e-2] <- 1e-2
        preds[preds>1-1e-2] <- 1-1e-2
        dataT <- cbind(dataT,preds)
      }
      data <- rbind(data,dataT)
    } 
    return(data)
  }
  fitAll <- dataCV[[1]]
  dataCV <- foreach (fold = unique(folds), .combine=cbind .inorder=FALSE) %dopar%
  {
    return(dataCV[[fold]])
  }
  
  #rename colnames of different lambda predictions
  NumberColumns <- 5
  colnames_lam <- paste0("lam1_",1:length(lam1V))
  colnames(dataCV)[1:length(lam1V) + NumberColumns] <- colnames_lam
  cv.results <- list(dataCV=dataCV)
  logLike <- matrix(0,nrow=length(tV),ncol=length(lam1V))
  colnames(logLike) <- colnames_lam
  rownames(logLike) <- paste0("TimePoint_",1:length(tV))
  AUC <- matrix(0,nrow=length(tV),ncol=length(lam1V))
  colnames(AUC) <- colnames_lam
  rownames(AUC) <- paste0("TimePoint_",1:length(tV))
  pWilcoxonMinusLog10 <- matrix(0,nrow=length(tV),ncol=length(lam1V))
  colnames(pWilcoxonMinusLog10) <- colnames_lam
  rownames(pWilcoxonMinusLog10) <- paste0("TimePoint_",1:length(tV))
  for (it in 1:length(tV))
  {
    IndT <- which(dataCV$timepoint==tV[it] & dataCV$status %in% c(0,1))
    predsT <- dataCV[IndT,colnames_lam]
    yT <- dataCV[IndT,"status"]
    weightsT <- (yT==0)/sum(yT==0) + (yT==1)/sum(yT==1)
    weightsT <- weightsT/sum(weightsT)
    for (ilam1 in 1:length(lam1V))
    {
      logLike[it,ilam1] <- sum(weightsT*yT*log(predsT[,ilam1]) + weightsT*(1-yT)*log(1-predsT[,ilam1]))
      AUC[it,ilam1] <- auc(yT, predsT[,ilam1], direction="<")[1]
      pWilcoxonMinusLog10[it,ilam1] <- -log10(wilcox.test(predsT[yT==1,ilam1], y = predsT[yT==0,ilam1], alternative = "greater", paired = FALSE, conf.int = FALSE)$p.value)
    }
  }
  if (whatToMaximize=="auc")
  {
    FigureMerit <- colMeans(AUC)
  }
  else if (whatToMaximize=="loglikelihood")
  {
    FigureMerit <- colMeans(logLike)
  }
  IndOptimum <- which.max(FigureMerit)
  OptimumLam1 <- lam1V[IndOptimum]
  return(list(dataCV=dataCV, logLike=logLike, AUC=AUC, pWilcoxonMinusLog10=pWilcoxonMinusLog10, OptimumLam1=lam1V[IndOptimum], gamma=gamma, fit=fitAll))
}





#' Fit lasso with non-penalized differences between adjacent time points coefficients using glmnet
#'
#' @param y0 case/control vector (no time iformation - "naive" approach)
#' @param x0 gene expression matrix (rows-samples by columns-genes)
#' @param FollowUp follow-up times (recurrence time for recurrences and follow-up for patients with no recurrences)
#' @param lam1V array of lasso penalty prefactor
#' @param gamma prefactor of the second penalty term - differences between adjacent time points coefficients
#' @param tV array of time points
#' @param Clinical0 dataframe with clinical information (same order as rows of x0)
#' @export
fitTimePointsNonPenalized <- function(y0, x0, FollowUp, lam1V, gamma, 
                                      tV, Clinical0=data.frame(case_control0=y0)){     
  y <- c()
  Clinical <- data.frame()
  samplesT <- 1:(nrow(x0)*length(tV))
  GenesT <- rep("",ncol(x0)*length(tV))
  Clinical0$sample <- rownames(x0)
  Clinical0$FollowUp <- FollowUp
  
  for (it in 1:length(tV)){
    t <- tV[it]
    case_controlT <- ifelse(FollowUp>t,0,ifelse(y0==1,1,-1))
    y <- c(y,case_controlT)
    ClinicalT <- Clinical0
    ClinicalT$StatusT <- case_controlT
    ClinicalT$FollowUp <- FollowUp
    ClinicalT$time <- tV[it]
    ClinicalT$samples <- rownames(x0)
    ClinicalT$samplesT <- paste0(rownames(x0),"_t_",it)
    Clinical <- rbind(Clinical,ClinicalT)
    GenesT[(1+(it-1)*ncol(x0)):(it*ncol(x0))] <- paste0(colnames(x0),"_t_",it)
    samplesT[(1+(it-1)*nrow(x0)):(it*nrow(x0))] <- paste0(rownames(x0),"_t_",it)
  }
  names(y) <- samplesT
  Ind <- which(y %in% c(0,1))
  Clinical <- Clinical[Ind,]
  y <- y[Ind]
  
  w <- y*0
  fits <- list()
  for (it in 1:length(tV)){
    IndT <- which(Clinical$time==tV[it])
    w[IndT][which(y[IndT]==0)] <- 1/sum(y[IndT]==0)
    w[IndT][which(y[IndT]==1)] <- 1/sum(y[IndT]==1)
    w[IndT] <- w[IndT]/sum(w)
    fit <- glmnet(x0[rownames(x0) %in% Clinical$samples[IndT],], y[IndT],
                  lambda=lam1V, intercept = TRUE, weights=w[IndT],  family = "binomial", alpha=1)
    fits <- list.append(fits,fit)
  }
  return(fits)
}       

