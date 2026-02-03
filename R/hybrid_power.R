#' Power Estimation for the multivariate goodness-of-fit problem via twosample tests
#' 
#' This function estimates the power of goodness-of-fit/two-sample hybrid tests using Rcpp and parallel computing
#' by generating a Monte Carlo data set and then running a twosample test.
#' 
#' For details on the usage of this routine consult the vignette with vignette("MDgof-hybrid","MDgof-hybrid")
#'
#' @param  rnull  routine to generate data under the null hypothesis.
#' @param  ralt  routine to generate data under the alternative hypothesis.
#' @param  param_alt values passed to ralt.
#' @param  phat =function(x) -99 parameter estimation, if needed.
#' @param  TS user supplied function to find test statistics, if any.
#' @param  TSextra (optional) list passed to TS, if needed.
#' @param  With.p.value =FALSE, does user supplied method find its own p-values?
#' @param  alpha =0.05 type I error rate used in tests.
#' @param  nMC =1 sample size of Monte Carlo data set, if it is a number nMC<=10  sample size used will be nMC*sample size of x.
#' @param  maxProcessor number of processors to use in parallel processing. 
#' @param  doMethods ="all", a vector of codes for the methods to include or all of them. 
#' @param  B =5000  number of simulation runs. If B=0 the routine returns the test statistics.
#' @return A list with vectors of test statistics and p.values
#' @examples
#' # All examples are run with B=20 and maxProcessor=1 to pass CRAN checks.
#' # Power of tests see whether data comes from a bivariate standard normal distribution, 
#' # without parameter estimation. True Distribution is bivariate normal with
#' # correlation r.
#' rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
#' ralt=function(r) mvtnorm::rmvnorm(100, sigma=matrix(c(1,r,r,1),2,2))
#' hybrid_power(rnull, ralt, 0.3, B=20, maxProcessor = 1)
#' # Power of tests to see whether data comes from a standard normal distribution, 
#' # with mean parameter estimated. True data comes from t distribution.
#' rnull=function(p) mvtnorm::rmvnorm(100, p)
#' ralt=function(df) mvtnorm::rmvt(100, df=df)
#' phat=function(x) apply(x, 2, mean)  
#' hybrid_power(rnull, ralt, 5, phat, B=20, maxProcessor = 1)
#' @export
hybrid_power=function(rnull, ralt, param_alt, phat=function(x) -99, 
                    nMC=1, TS, TSextra, With.p.value=FALSE,
                    alpha=0.05, B=1000,
                    maxProcessor, doMethods="all") {
  x=ralt(param_alt[1]) # example data set
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
  f=function(a) {
    if(length(formals(rnull))==0)  {
      xn=rnull()
      yn=ralt(a)
    }
    else {
      xn=rnull(p)
      yn=ralt(a)
    }
    if(Continuous) {
      while (nrow(yn)<nMC) yn=rbind(yn, ralt(a))
      yn=yn[1:nMC, ]
      out=list(x=xn, y=yn)
    }
    else {
      while (sum(yn[,3])<nMC) yn[,3]=yn[,3]+ralt(a)[,3]
      out=cbind(xn[,1:3], yn[,3])
    }   
    out
  }
  MD2sample::twosample_power(f, param_alt, 
        TS=TS, TSextra=TSextra, alpha=alpha, 
        maxProcessor=maxProcessor, 
        With.p.value=With.p.value, doMethods=doMethods, B=B)
}
