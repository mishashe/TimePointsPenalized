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
#' @param standardize TRUE/FALSE standardization of the x0 columns (zero mean, unit variance)
#' @param Clinical0 dataframe with clinical information (same order as rows of x0)
#' @export
fitTimePointsPenalized <- function(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinical0=data.frame(case_control0=y0))
{     
  if (standardize) {
    for (i in 1:ncol(x0)) {
      x0[,i] <- (x0[,i] - mean(x0[,i]))/sd(x0[,i])
    }
  }
  Intercept <- rep(0,length(tV))
  beta <- rep(0,ncol(x0)*length(tV))
  y <- c()
  Clinical <- data.frame()
  samplesT <- 1:(nrow(x0)*length(tV))
  GenesT <- rep("",ncol(x0)*length(tV))
  Clinical0$sample <- rownames(x0)
  Clinical0$FollowUp <- FollowUp
  
  for (it in 1:length(tV))
  {
    t <- tV[it]
    case_controlT <- ifelse(FollowUp>t,0,ifelse(y0==1,1,-1))
    if (sum(case_controlT==0)<2)
    {
      print(paste0("Not enough controls for t=",t))
      print(table(case_controlT))
      return(NULL)
    }
    if (sum(case_controlT==1)<2)
    {
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
  for (it in 1:length(tV))
  {
    IndT <- which(Clinical$time==tV[it])
    w[IndT][which(y[IndT]==0)] <- 1/sum(y[IndT]==0)
    w[IndT][which(y[IndT]==1)] <- 1/sum(y[IndT]==1)
    IndFor0 <- c(IndFor0,which(rownames(x0) %in% Clinical$samples[IndT]))
    IndTFor0 <- c(IndTFor0,which(rownames(x0) %in% Clinical$samples[IndT])*0+it)
  }
  
  fits <- list(list())
  for (it in 1:length(tV))
  {
    fits[[it]] <- list()
    fits[[it]]$beta <- matrix(0,nrow=dim(x0)[2],ncol=0)
    fits[[it]]$Intercept <- c()
  }
  for (ilam1 in 1:length(lam1V))
  {
    lam1 <- lam1V[ilam1]
    lam2 <- gamma*lam1
    
    fit <- Fit(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0, IndTFor0)
    
    for (it in 1:length(tV))
    {
      print(fit$beta)
      beta <- fit$beta[(1:(dim(x0)[2]))+(it-1)*dim(x0)[2],1,drop=TRUE]
      print(2)
      Intercept <- fit$Intercept[it]
      print(3)
      IndT <- which(Clinical$time==tV[it])
      print(4)
      fits[[it]]$beta <- cbind(fits[[it]]$beta,beta)
      print(5)
      fits[[it]]$Intercept <- c(fits[[it]]$Intercept,Intercept)
      print(6)
    }
  }
  
  for (it in 1:length(tV))
  {
    IndT <- which(Clinical$time==tV[it])
    
  }
  return(fitsT)
}       

#' Fit lasso with non-penalized differences between adjacent time points coefficients using glmnet
#'
#' @param y0 case/control vector (no time iformation - "naive" approach)
#' @param x0 gene expression matrix (rows-samples by columns-genes)
#' @param FollowUp follow-up times (recurrence time for recurrences and follow-up for patients with no recurrences)
#' @param lam1V array of lasso penalty prefactor
#' @param gamma prefactor of the second penalty term - differences between adjacent time points coefficients
#' @param tV array of time points
#' @param standardize TRUE/FALSE standardization of the x0 columns (zero mean, unit variance)
#' @param Clinical0 dataframe with clinical information (same order as rows of x0)
#' @export
fitTimePointsNonPenalized <- function(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinical0=data.frame(case_control0=y0))
{     
  if (standardize) {
    for (i in 1:ncol(x0)) {
      x0[,i] <- (x0[,i] - mean(x0[,i]))/sd(x0[,i])
    }
  }

  y <- c()
  Clinical <- data.frame()
  samplesT <- 1:(nrow(x0)*length(tV))
  GenesT <- rep("",ncol(x0)*length(tV))
  Clinical0$sample <- rownames(x0)
  Clinical0$FollowUp <- FollowUp
  
  for (it in 1:length(tV))
  {
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
  names(beta) <- GenesT
  names(y) <- samplesT
  Ind <- which(y %in% c(0,1))
  Clinical <- Clinical[Ind,]
  y <- y[Ind]
  
  IndFor0 <- c()
  IndTFor0 <- c()
  w <- y*0
  for (it in 1:length(tV))
  {
    IndT <- which(Clinical$time==tV[it])
    w[IndT][which(y[IndT]==0)] <- 1/sum(y[IndT]==0)
    w[IndT][which(y[IndT]==1)] <- 1/sum(y[IndT]==1)
    IndFor0 <- c(IndFor0,which(rownames(x0) %in% Clinical$samples[IndT]))
    IndTFor0 <- c(IndTFor0,which(rownames(x0) %in% Clinical$samples[IndT])*0+it)
  }
  fits <- list()
  
  for (it in 1:length(tV))
  {
    IndT <- which(Clinical$time==tV[it])
    w[IndT][which(y[IndT]==0)] <- 1/sum(y[IndT]==0)
    w[IndT][which(y[IndT]==1)] <- 1/sum(y[IndT]==1)
    fit <- glmnet(x0[which(rownames(x0) %in% Clinical$samples[IndT]),], y[IndT],
                  lambda=lam1V, standardize = FALSE, intercept = TRUE, weights=w[IndT],  family = "binomial", alpha=1)
    fits <- list.append(fits,fit)
  }
  return(fits)
}       
