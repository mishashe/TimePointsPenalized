library(rlist) library(Rcpp) library(devtools)
detach(“package:TimePointsPenalized”, unload=TRUE)
remove.packages(“TimePointsPenalized”)
install\_github(“mishashe/TimePointsPenalized”, force=TRUE)
library(TimePointsPenalized) nSamples &lt;- 100 nGenes &lt;- 300 tV
&lt;- seq(2,8,1)*12 x0 &lt;-
matrix(rnorm(nSamples*nGenes),nrow=nSamples,ncol=nGenes) beta0 &lt;-
rep(0,nGenes) beta0\[1:5\] &lt;- 1 beta0\[6:10\] &lt;- -1 y0 &lt;-
round(1/(1+exp(-x0 %*% beta0))) beta &lt;- rep(0,nGenes*length(tV)) tV
&lt;- seq(2,8,1)*12 FollowUp &lt;- 2*12 + (8-2)*12*runif(nSamples) lam1V
&lt;- 10^seq(-1.0,-4.5,-0.025) gamma &lt;- 10 rownames(x0) &lt;-
paste0(“S”,1:nSamples) colnames(x0) &lt;- paste0(“G”,1:nGenes) fits
&lt;- fitTimePointsPenalized(y0, x0, FollowUp, lam1V, gamma, tV,
standardize=TRUE, Clinical0=data.frame(case\_control0=y0))
