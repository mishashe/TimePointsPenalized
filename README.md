```
library(rlist)
library(Rcpp)
library(devtools)
detach("package:TimePointsPenalized", unload=TRUE)
remove.packages("TimePointsPenalized")
install_github("mishashe/TimePointsPenalized", force=TRUE)
library(TimePointsPenalized)
nSamples <- 100
nGenes <- 300
tV <- seq(2,8,1)*12
x0 <- matrix(rnorm(nSamples*nGenes),nrow=nSamples,ncol=nGenes)
beta0 <- rep(0,nGenes)
beta0[1:5] <- 1
beta0[6:10] <- -1
y0 <- round(1/(1+exp(-x0 %*% beta0)))
beta <- rep(0,nGenes*length(tV))
tV <- seq(2,8,1)*12
FollowUp <- 2*12 + (8-2)*12*runif(nSamples)
lam1V <- 10^seq(-1.0,-4.5,-0.025)
gamma <- 10
rownames(x0) <- paste0("S",1:nSamples)
colnames(x0) <- paste0("G",1:nGenes)
fits <- fitTimePointsPenalized(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinical0=data.frame(case_control0=y0))
```
