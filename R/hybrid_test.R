#' Tests for the multivariate goodness-of-fit problem via twosample tests
#' 
#' This function runs a number of goodness-of-fit tests using Rcpp and parallel computing
#' by generating a Monte Carlo data set and then running a twosample test.
#' 
#' For details on the usage of this routine consult the vignette with vignette("MDgof-hybrid","MDgof-hybrid")
#'
#' @param  x a matrix with the data set
#' @param  rnull  routine to generate data under the null hypothesis.
#' @param  phat =function(x) -99 parameter estimation, if needed.
#' @param  TS user supplied function to find test statistics, if any.
#' @param  TSextra (optional) list passed to TS, if needed.
#' @param  nMC =1 sample size of Monte Carlo data set, if it is a number nMC<=10  sample size used will be nMC*sample size of x.
#' @param  maxProcessor number of processors to use in parallel processing. 
#' @param  doMethods ="all", a vector of codes for the methods to include or all of them. 
#' @param  B =5000  number of simulation runs. If B=0 the routine returns the test statistics.
#' @return A list with vectors of test statistics and p.values
#' @examples
#' # All examples are run with B=20 and maxProcessor=1 to pass CRAN checks.
#' # Tests to see whether data comes from a bivariate standard normal distribution, 
#' # without parameter estimation.
#' rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
#' x=rnull()
#' hybrid_test(x, rnull, B=20, maxProcessor = 1)
#' # Tests to see whether data comes from a standard normal distribution, 
#' # with mean parameter estimated.
#' rnull=function(p) mvtnorm::rmvnorm(100, p)
#' phat=function(x) apply(x, 2, mean)  
#' x=rnull(c(0,1))
#' hybrid_test(x, rnull, phat, B=20, maxProcessor = 1)
#' # Example of a discrete model, without parameter estimation
#' # X~Bin(5, 0.5), Y|X=x~Bin(4, 0.5+x/100)
#' rnull=function() {
#'   x=rbinom(1000, 5, 0.5)
#'   y=rbinom(1000, 4, 0.5)
#'   MDgof::sq2rec(table(x, y))
#' }
#' x=rnull()
#' hybrid_test(x, rnull, B=50, maxProcessor = 1)
#' # Example of a discrete model, with parameter estimation
#' # X~Bin(5, p), Y|X=x~Bin(4, 0.5+x/100)
#' rnull=function(p) {
#'   x=rbinom(1000, 5, p)
#'   y=rbinom(1000, 4, 0.5+x/100)
#'   MDgof::sq2rec(table(x, y))
#' }
#' phat=function(x) {
#'   tx=tapply(x[,3], x[,1], sum)
#'   p1=mean(rep(as.numeric(names(tx)), times=tx))/5
#'   p1
#' }
#' x=rnull(0.5)
#' hybrid_test(x, rnull, phat, B=20, maxProcessor = 1)
#' @export
hybrid_test=function(x, rnull, phat=function(x) -99, nMC=1, 
                     TS, TSextra,  
                     B=1000, maxProcessor, doMethods="all") {
  Continuous=TRUE
  if(ncol(x)==3 && all(abs(round(x[,3])-x[,3])<0.001))
      Continuous=FALSE
  if(Continuous) {
      p=phat(x)
      if(nMC<=10) nMC=nMC*nrow(x)
  }
  else {
      p=phat(x)
      if(nMC<=10) nMC=nMC*sum(x[,3])
  }
  nMC=round(nMC)
  rnullT=function(dta) {
    if(length(formals(rnull))==0)  {
      xn=rnull()
      yn=rnull()
    }
    else {
      xn=rnull(p)
      yn=rnull(p)
    }
    if(Continuous) {
      while (nrow(yn)<nMC) {
        if(length(formals(rnull))==0)  yn=rbind(yn, rnull())
        else yn=rbind(yn, rnull(p))
      }
      yn=yn[1:nMC, ]
      out=list(x=xn, y=yn)
    }
    else {
      while (sum(yn[,3])<nMC) {
        if(length(formals(rnull))==0)  yn[,3]=yn[,3]+rnull()[,3]
        else yn[,3]=yn[,3]+rnull(p)[,3]
      }
      out=cbind(xn[,1:3], yn[,3])
      colnames(out)=c("vals_x", "vals_y", "x", "y")
    }   
    out
  }
  if(Continuous)
    out=MD2sample::twosample_test(x, rnullT(list(x=x))$y, 
              TS=TS, TSextra=TSextra, 
              rnull=rnullT, maxProcessor=maxProcessor, 
              doMethods=doMethods, B=B)
  else 
    out=MD2sample::twosample_test(x[,3], rnullT(x)[,4], 
              x[,1], x[,2], TS=TS, TSextra=TSextra, 
              rnull=rnullT, maxProcessor=maxProcessor, 
              doMethods=doMethods, B=B)
  out
}
