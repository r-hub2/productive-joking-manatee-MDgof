#' Power estimation of goodness-of-fit tests.
#' 
#' Find the power of various goodness-of-fit tests.
#' 
#' For details on the usage of this routine consult the vignette with vignette("MDgof","MDgof")
#' 
#' @param  pnull function to find cdf under  null hypothesis
#' @param  rnull function to generate data under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99
#' @param  dnull =function(x) -99, density function under the null hypothesis, if available, or -99 if missing
#' @param  TS user supplied function to find test statistics
#' @param  TSextra list provided to TS (optional)
#' @param  With.p.value =FALSE does user supplied routine return p values?
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2), a 2x2 matrix with lower and upper bounds, if any, for chi-square tests
#' @param  nbins =c(5, 5), number of bins for chi square tests.
#' @param  minexpcount =5 minimal expected bin count required for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  SuppressMessages =FALSE, should informative messages be printed?
#' @param  maxProcessor maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  B =1000 number of simulation runs
#' @return A numeric matrix of power values.
#' @examples
#' # All examples are run with B=10 and maxProcessor=1 to pass CRAN checks.
#' # This is obviously MUCH TO SMALL for any real usage.
#' # Power of tests if null hypothesis specifies a bivariate standard normal 
#' # distribution but data comes from a bivariate normal with different means, 
#' # without parameter estimation.
#' rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
#' ralt=function(p) mvtnorm::rmvnorm(100, c(p, p))
#' pnull=function(x) {
#'   if(!is.matrix(x)) return(mvtnorm::pmvnorm(rep(-Inf, 2), x))
#'   apply(x, 1, function(x) mvtnorm::pmvnorm(rep(-Inf, 2), x))
#' }
#' gof_power(pnull, rnull, ralt, c(0, 1), B=10, maxProcessor = 1)
#' # Same as above, but now with density included
#' dnull=function(x) {
#'   if(!is.matrix(x)) return(mvtnorm::dmvnorm(x))
#'   apply(x, 1, function(x) mvtnorm::dmvnorm(x))
#' }
#' \donttest{gof_power(pnull, rnull, ralt, c(0, 1), dnull=dnull, B=10, maxProcessor = 1)}
#' # Power of tests when null hypothesis specifies a bivariate normal distribution, 
#' # with mean parameter estimated, wheras data comes from a t distribution
#' rnull=function(p) mvtnorm::rmvnorm(100, p)
#' ralt=function(df) mvtnorm::rmvt(100, sigma=diag(2), df=df)
#' pnull=function(x,p) {
#'   if(!is.matrix(x)) return(mvtnorm::pmvnorm(rep(-Inf, 2), x, mean=p))
#'   apply(x, 1, function(x) mvtnorm::pmvnorm(rep(-Inf, 2), x, mean=p))
#' }
#' dnull=function(x, p) {
#'   if(!is.matrix(x)) return(mvtnorm::dmvnorm(x, mean=p))
#'   apply(x, 1, function(x) mvtnorm::dmvnorm(x, mean=p))
#' }
#' phat=function(x) apply(x, 2, mean)
#' \donttest{gof_power(pnull, rnull, ralt, c(50, 5), dnull=dnull, phat=phat, B=10, maxProcessor = 1)}
#' # Example of a discrete model, with parameter estimation
#' # Under null hypothesis: X~Bin(10, p), Y|X=x~Bin(5, 0.5+x/100)
#' # Under alternative hypothesis: X~Bin(10, p), Y|X=x~Bin(5, K+x/100)
#' rnull=function(p=0.5) {
#'   x=stats::rbinom(1000, 10, p)
#'   y=stats::rbinom(1000, 5, 0.5+x/100)
#'   MDgof::sq2rec(table(x, y))
#' }
#' ralt=function(K=0.5) {
#'   x=stats::rbinom(1000, 10, 0.5)
#'   y=stats::rbinom(1000, 5, K+x/100)
#'   MDgof::sq2rec(table(x, y))
#' }
#' pnull=function(x, p) {
#'   f=function(x) sum(dbinom(0:x[1], 10, p[1])*pbinom(x[2], 5, 0.5+0:x[1]/100))
#'   if(!is.matrix(x)) x=rbind(x)
#'   apply(x, 1, f)
#' }
#' phat=function(x) {
#'   tx=tapply(x[,3], x[,1], sum)
#'   mean(rep(as.numeric(names(tx)), times=tx))/10
#' }
#' \donttest{gof_power(pnull, rnull, ralt, c(0.5, 0.6), phat=phat, B=10, maxProcessor = 1)}
#' @export 
gof_power=function(pnull, rnull, ralt, param_alt, 
        phat=function(x) -99, dnull=function(x) -99, TS, TSextra, 
        With.p.value=FALSE, 
        alpha=0.05, Ranges=matrix(c(-Inf, Inf, -Inf, Inf),2,2), nbins=c(5 ,5), 
        minexpcount=5.0, rate=0, SuppressMessages=FALSE,maxProcessor,  B=1000) {
  NewTS=ifelse(missing(TS), FALSE, TRUE)
  dta = ralt(param_alt[1]) # get an example data set
  Continuous=TRUE
  if(ncol(dta)==3 && (all(abs(round(dta[,3])-dta[,3])<0.001))) 
     Continuous=FALSE
  check.functions(pnull, rnull, phat, x=dta) 
  TSextra=makeTSextra(dta, Continuous, pnull, rnull, phat, dnull, Ranges, TSextra)
  #If new test finds p values, run power study now
  if(With.p.value) {
    out=power_pvals(pnull, ralt, param_alt, TS, TSextra, alpha, B)
    return(round(out, 3))
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
      if(!SuppressMessages)
        message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessor=1
    }
    if(length(formals(TS))<3 | length(formals(TS))>4) {
      message("TS should have either 3 or 4 arguments (dta, pnull, param and TSextra (optional)")
      return(NULL)
    }
    typeTS=ifelse(length(formals(TS))==3, 1, 2)
  }
  TS_data=calcTS(dta, TS, typeTS, TSextra) 
  if(is.null(names(TS_data))) {
    message("result of TS has to be a named vector")
    return(NULL)
  }
  if(missing(maxProcessor)) 
    maxProcessor=parallel::detectCores(logical = FALSE)-1
    if(With.p.value) maxProcessor=1
  if(maxProcessor>1) {
    tm=timecheck(dta, TS, typeTS, TSextra)
    if(tm*length(param_alt)*B<20) {
      maxProcessor=1
      if(!SuppressMessages)
          message("maxProcessor set to 1 for faster computation")
    }
    else {
      if(!SuppressMessages)  
         message(paste("Using ",maxProcessor," cores.."))
    }  
  }
  if(maxProcessor==1) {
    tmp=powerC(rnull, ralt, param_alt, TS, typeTS, TSextra, B)
    Data=tmp$Data
    Sim=tmp$Sim
    param_alt_vec=tmp$param_alt_vec
  }
  else {
    cl <- parallel::makeCluster(maxProcessor)
    z=parallel::clusterCall(cl, powerC, 
                rnull, ralt, param_alt,  TS, typeTS, TSextra, 
                B=round(B/maxProcessor))
    parallel::stopCluster(cl)
    Sim=z[[1]][["Sim"]]
    Data=z[[1]][["Data"]]
    param_alt_vec=z[[1]][["param_alt_vec"]]
    for(i in 2:maxProcessor) {
        Sim=rbind(Sim,z[[i]][["Sim"]])
        Data=rbind(Data,z[[i]][["Data"]])
        param_alt_vec=c(param_alt_vec, z[[i]][["param_alt_vec"]])
    }  
  }  
  pwr=matrix(0, length(param_alt), length(TS_data))
  colnames(pwr)=names(TS_data)
  rownames(pwr)=param_alt
  crtval=apply(Data, 2, quantile, prob=1-alpha, na.rm=TRUE)
  for(i in seq_along(param_alt)) {
    tmpS=Sim[param_alt_vec==param_alt[i], ,drop=FALSE]
    for(j in seq_along(crtval)) 
       pwr[i, j]=sum(tmpS[ ,j]>crtval[j])/nrow(tmpS)
  }
  # Do chi square tests if built-in TS is used, data is continuous and Dim=2. 
  chipwr=NULL
  if(!NewTS & ncol(dta)==2) {
      chipwr = chi_power(pnull=pnull, 
                         ralt = ralt, 
                         param_alt = param_alt,                             
                         phat = phat, 
                         alpha = alpha, 
                         Ranges = Ranges, 
                         nbins = nbins, 
                         rate = rate, 
                         minexpcount = minexpcount,
                         B = B) 
  }
  pwr = cbind(pwr, chipwr)
  if(TSextra$NoDensity) {
    pwr=pwr[ ,colnames(pwr)!="BB", drop=FALSE]
    pwr=pwr[ ,colnames(pwr)!="KSD", drop=FALSE]
  }
  if(ncol(dta)>2) {
    pwr=pwr[ ,colnames(pwr)!="MMD", drop=FALSE]
    pwr=pwr[ ,colnames(pwr)!="FF", drop=FALSE]
    pwr=pwr[ ,colnames(pwr)!="RK", drop=FALSE]
  }  
  if(is.matrix(pwr) & nrow(pwr)==1) pwr=pwr[1, ]
  round(pwr, 3)
}
