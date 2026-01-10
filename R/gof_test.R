#' Tests for the multivariate goodness-of-fit problem
#' 
#' This function runs a number of goodness-of-fit tests using Rcpp and parallel computing.
#' 
#' For details on the usage of this routine consult the vignette with vignette("MDgof","MDgof")
#'
#' @param  x a matrix with the data set
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  phat  =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  dnull =function(x) -99, density function under the null hypothesis, if available, or -99 if missing
#' @param  TS user supplied function to find test statistics, if any.
#' @param  TSextra (optional) list passed to TS, if needed.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  nbins =c(5, 5) number of bins for chi-square tests
#' @param  Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2), a 2x2 matrix with lower and upper bounds, if any, for chi-square tests
#' @param  minexpcount =5 minimal expected bin count required
#' @param  maxProcessor number of processors to use in parallel processing. 
#' @param  doMethods ="all", a vector of codes for the methods to include or all of them. 
#' @param  B =5000  number of simulation runs. If B=0 the routine returns the test statistics.
#' @param  ReturnTSextra =FALSE, should setup info be returned?
#' @return A list with vectors of test statistics and p.values
#' @examples
#' # All examples are run with B=20 and maxProcessor=1 to pass CRAN checks.
#' # Tests to see whether data comes from a bivariate standard normal distribution, 
#' # without parameter estimation.
#' rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
#' x=rnull()
#' pnull=function(x) {
#'   if(!is.matrix(x)) return(mvtnorm::pmvnorm(rep(-Inf, 2), x))
#'   apply(x, 1, function(x) mvtnorm::pmvnorm(rep(-Inf, 2), x))
#' }
#' gof_test(x, pnull, rnull, B=20, maxProcessor = 1)
#' # Same as above, but now with density included
#' dnull=function(x) {
#'   if(!is.matrix(x)) return(mvtnorm::dmvnorm(x))
#'   apply(x, 1, function(x) mvtnorm::dmvnorm(x))
#' }
#' gof_test(x, pnull, rnull, dnull=dnull, B=20, maxProcessor = 1)
#' # Tests to see whether data comes from a standard normal distribution, 
#' # with mean parameter estimated.
#' rnull=function(p) mvtnorm::rmvnorm(100, p)
#' x=rnull(c(0,1))
#' pnull=function(x,p) {
#'   if(!is.matrix(x)) return(mvtnorm::pmvnorm(rep(-Inf, 2), x, mean=p))
#'   apply(x, 1, function(x) mvtnorm::pmvnorm(rep(-Inf, 2), x, mean=p))
#' }
#' dnull=function(x, p) {
#'   if(!is.matrix(x)) return(mvtnorm::dmvnorm(x, mean=p))
#'   apply(x, 1, function(x) mvtnorm::dmvnorm(x, mean=p))
#' }
#' phat=function(x) apply(x, 2, mean)
#' gof_test(x, pnull, rnull, dnull=dnull, phat=phat,B=20, maxProcessor = 1)
#' # Example of a discrete model, with parameter estimation
#' # X~Bin(10, p1), Y|X=x~Bin(5, p2+x/100)
#' rnull=function(p) {
#'   x=rbinom(1000, 10, p[1])
#'   y=rbinom(1000, 5, p[2]+x/100)
#'   MDgof::sq2rec(table(x, y))
#' }
#' pnull=function(x, p) {
#'   f=function(x) sum(dbinom(0:x[1], 10, p[1])*pbinom(x[2], 5, p[2]+0:x[1]/100))
#'   if(!is.matrix(x)) x=rbind(x)
#'   apply(x, 1, f)
#' }
#' phat=function(x) {
#'   tx=tapply(x[,3], x[,1], sum)
#'   p1=mean(rep(as.numeric(names(tx)), times=tx))/10
#'   ty=tapply(x[,3], x[,2], sum)
#'   p2=mean(rep(as.numeric(names(ty)), times=ty))/5-p1/10
#'   c(p1, p2)
#' }
#' x=rnull(c(0.5, 0.5))
#' gof_test(x, pnull, rnull, phat=phat,B=20, maxProcessor = 1)
#' @export
gof_test <- function(x, pnull, rnull, phat=function(x) -99, 
                    dnull=function(x) -99,
                    TS, TSextra, rate=0,  nbins=c(5, 5),  
                    Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2), minexpcount=5.0,
                    maxProcessor, doMethods="all", B=5000,
                    ReturnTSextra=FALSE) {
  check.functions(pnull, rnull, phat, x=x)
  NewTS=ifelse(missing(TS), FALSE, TRUE)
  Continuous=TRUE
  if(ncol(x)==3 && all(abs(round(x[,3])-x[,3])<0.001)) Continuous=FALSE
  if(doMethods[1]=="all") {
    if(Continuous) {
       doMethods=c("qKS", "qK", "qCvM", "qAD", "BB", "BR",  "KSD")
       if(ncol(x)==2)
          doMethods=c(doMethods, "MMD", "FF", "RK", "ES", "EP")
    }   
    if(!Continuous)
        doMethods=c("KS", "K", "CvM", "AD", "TV", "KL", "H", "P")
  }
  TSextra=makeTSextra(x, Continuous, pnull, rnull, phat, dnull, Ranges, TSextra) 
  if(ReturnTSextra) return(TSextra)
  if(TSextra$NoDensity) {
    doMethods=doMethods[doMethods!="BB"]
    doMethods=doMethods[doMethods!="KSD"]
  }
# Set up some things:
   if(!NewTS) {
      typeTS=2
      if(Continuous) TS = TS_cont
      else TS = TS_disc
   }   
   else {
# can't do parallel processing if TS written in C/C++
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
       message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
       maxProcessor=1
    }
    if(length(formals(TS))<3 | length(formals(TS))>4) {
      message("TS should have either 3 or 4 arguments (dta, pnull, param and TSextra (optional)")
      return(NULL)
    }
    typeTS=ifelse(length(formals(TS))==3, 1, 2)
  }
  TS_data=calcTS(x, TS, typeTS, TSextra)  
  if(B==0) return(TS_data)
  if(is.null(names(TS_data))) {
    message("result of TS has to be a named vector")
    return(NULL)
  }
  if(missing(maxProcessor)) 
     maxProcessor=parallel::detectCores(logical = FALSE)-1
  if(maxProcessor>1) {
     tm=timecheck(x, TS, typeTS, TSextra)
     if(tm*B<20) {
       maxProcessor=1
       message("maxProcessor set to 1 for faster computation")
     }
     else message(paste("Using ",maxProcessor," cores.."))
  }
  if(maxProcessor==1) outTS=testC(x, rnull, TS, typeTS, TSextra, B) 
  else {
     cl = parallel::makeCluster(maxProcessor)
     z=parallel::clusterCall(cl, testC, 
             x, rnull, TS, typeTS, TSextra, 
             B=round(B/maxProcessor)
     )
     parallel::stopCluster(cl)
     #  Average power of cores
     tmp=0*z[[1]]
     for(i in 1:maxProcessor) tmp=tmp+z[[i]]
     outTS = tmp/maxProcessor  
  }
# do chi square tests if Continuous and Dim=2
  outchi=NULL
  if(!NewTS & typeTS==2) {
     if(Continuous & ncol(x)==2) {
       tmp=chi_cont_test(x, pnull, phat, Ranges, nbins, minexpcount)
       outchi = t(tmp[,c(1, 2)])
     }
  }
  out=list(statistics=c(outTS[1, ], outchi[1, ]), 
             p.values=c(outTS[2, ], outchi[2, ]))
  if(!NewTS) {
    out[[1]]=out[[1]][doMethods]
    out[[2]]=out[[2]][doMethods]
  }  
  # make output look nice
  signif.digits(out)
}
