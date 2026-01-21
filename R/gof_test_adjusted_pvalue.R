#' Helper function to find test statistics of simulated data.
#' @param dta a matrix with data
#' @param TS test statistic routine
#' @param typeTS type of routine
#' @param TSextra a list
#' @param B number of simulation runs
#' @return a matrix
#' @keywords internal
#' @export
simTS=function(dta, TS, typeTS, TSextra, B) {
  A=matrix(0, B, length(calcTS(dta, TS, typeTS, TSextra)))
  p=TSextra$phat(dta)
  for(i in 1:B) {
    if(length(formals(TSextra$rnull))==0)
      simdta=TSextra$rnull()
    else 
      simdta=TSextra$rnull(p)
    if(!TSextra$NoDensity & TSextra$Continuous) 
      TSextra$nn=FNN::get.knn(simdta, 1)
    A[i,]=calcTS(simdta, TS, typeTS, TSextra)
  }
  A
}

#' Helper function to find p values of simulated data.
#' @param  dta a matrix with data
#' @param TS test statistic routine
#' @param typeTS type of routine
#' @param TSextra a list
#' @param A a matrix
#' @param Ranges a matrix
#' @param nbins a vector
#' @param minexpcount an integer
#' @param B number of simulation runs
#' @return a matrix
#' @keywords internal
#' @export
simpvals=function(dta, TS, typeTS, TSextra, A,  
                  Ranges, nbins, minexpcount, B) {
  num_tests=length(calcTS(dta, TS, typeTS, TSextra))
  pvalsTS=matrix(0, B, num_tests)
  pvalsChi=NULL
  if(TSextra$Continuous) {
       pvalsChi=matrix(0, B, 2)
       colnames(pvalsChi)=c("ES","EP")
  }   
  p=TSextra$phat(dta)  
  for(i in 1:B) {
      if(length(formals(TSextra$rnull))==0)
        simdta=TSextra$rnull()
      else 
        simdta=TSextra$rnull(p)
    if(!TSextra$NoDensity & TSextra$Continuous) 
        TSextra$nn=FNN::get.knn(simdta, 1)  
    tmp=calcTS(simdta, TS, typeTS, TSextra)
    if(TSextra$Continuous) {
        pvalsChi[i, ]=chi_cont_test(simdta, TSextra$pnull, 
                                    TSextra$phat, Ranges, nbins, minexpcount)[,2]
    }
    for(j in 1:num_tests) pvalsTS[i,j]=pvalsTS[i,j]+sum(tmp[j]>A[,j])/nrow(A) 
  }
  colnames(pvalsTS)=names(tmp)
  list(pvalsTS=pvalsTS, pvalsChi=pvalsChi)
}

#' Adjusted p values
#' 
#' This function runs a number of goodness-f-fit tests using Rcpp and parallel computing and then finds the correct p value for the combined tests.
#' 
#' For details consult the vignette("MDgof","MDgof")
#' 
#' @param  x  matrix with data
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  phat  =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  dnull =function(x) -99, density function under the null hypothesis, if available, or -99 if missing
#' @param  B =c(5000, 1000), number of simulation runs for permutation test and for estimation
#'         of the empirical distribution function.
#' @param  nbins =c(5, 5), number of bins for chi square tests (2D only).
#' @param  minexpcount = 5, minimum required expected counts for chi-square tests.
#' @param  Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2) a 2x2 matrix with lower and upper bounds.
#' @param  SuppressMessages = FALSE, print informative messages?
#' @param  maxProcessor number of cores for parallel processing.
#' @param  doMethods  Which methods should be included? If missing a small number of methods that generally have good power are used.
#' @return NULL, results are printed out.
#' @examples
#' # All examples are run with B=10 and maxProcessor=1 to pass CRAN checks.
#' # This is obviously MUCH TO SMALL for any real usage.
#' # Tests to see whether data comes from a bivariate standard normal distribution, 
#' # without parameter estimation.
#' rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
#' x=rnull()
#' pnull=function(x) {
#'   if(!is.matrix(x)) return(mvtnorm::pmvnorm(rep(-Inf, 2), x))
#'   apply(x, 1, function(x) mvtnorm::pmvnorm(rep(-Inf, 2), x))
#' }
#' dnull=function(x) {
#'   if(!is.matrix(x)) return(mvtnorm::dmvnorm(x))
#'   apply(x, 1, function(x) mvtnorm::dmvnorm(x))
#' }
#' gof_test_adjusted_pvalue(x, pnull, rnull, dnull=dnull, B=10, maxProcessor = 1)
#' @export
gof_test_adjusted_pvalue=function(x, pnull, rnull, 
                        phat=function(x) -99,  
                        dnull=function(x) -99,  
                        B=c(5000, 1000), nbins=c(5,5),
                        minexpcount=5, 
                        Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2),
                        SuppressMessages=FALSE, 
                        maxProcessor, doMethods="all") {
  Continuous=TRUE
  if(ncol(x)==3 && all(abs(round(x[,3])-x[,3])<0.001)) Continuous=FALSE
  if(length(B)==1) B=c(B, B)
  TSextra=makeTSextra(x, Continuous, pnull, rnull, phat, dnull, Ranges)   
  if(Continuous) {
     allMethods=c("qKS", "qK", "qCvM", "qAD", "BB", "BR",  "KSD")
     if(ncol(x)==2)  allMethods=c(allMethods, "MMD", "FF", "RK", "ES", "EP")
     if(TSextra$NoDensity) 
       IndexNotIncluded=seq_along(allMethods)[allMethods=="BB" | allMethods=="KSD" ]
  }   
  if(!Continuous) allMethods=c("KS", "K", "CvM", "AD", "TV", "KL", "H", "P")
  if(doMethods[1]=="all") {
     doMethods=allMethods
     if(TSextra$NoDensity) doMethods=doMethods[-IndexNotIncluded]
  }   
  else {
    if(TSextra$NoDensity & "BB"%in%doMethods) {
       message("Method BB can not be run without dnull (density function)")
      doMethods=doMethods[doMethods!="BB"]
    }
    if(TSextra$NoDensity & "KSD"%in%doMethods) {
      message("Method KSD can not be run without dnull (density function)")
      doMethods=doMethods[doMethods!="KSD"]
    }
  }
  
  if(Continuous) {
    if(ncol(x)==2) {
      if(length(nbins)==1) nbins=c(nbins, nbins)
      outchi = MDgof::chi_cont_test(x, pnull, phat, Ranges, nbins, minexpcount)[,2]
    }  
    typeTS=2
    TS=TS_cont
  }
  else {
    typeTS=2
    outchi = MDgof::chi_disc_test(x, pnull, dnull, phat, minexpcount)[2]
    TS=TS_disc
  }    
  TS_data=calcTS(x, TS, typeTS, TSextra)
  if(any(is.null(names(TS_data)))) {
    if(!SuppressMessages) message("output of TS routine has to be a named vector!")
    return(NULL)
  }  
#   Find matrix of values of test statistics for simulated data
  if(missing(maxProcessor)) {
     maxProcessor=parallel::detectCores(logical = FALSE)-1
     if(!SuppressMessages) message(paste("Using ", maxProcessor," cores.."))
  }   
  if(maxProcessor==1) A=simTS(x, TS, typeTS, TSextra, B[1])
  else {
     cl=parallel::makeCluster(maxProcessor)
     z=parallel::clusterCall(cl, simTS, x, TS, typeTS, 
                               TSextra, round(B[1]/maxProcessor))
     A=z[[1]]
     for(i in 2:maxProcessor) A=rbind(A, z[[i]])
     B[1]=nrow(A)
  }
  num_tests=ncol(A)
  pvalsdta=rep(0, num_tests)
  for(j in 1:num_tests) pvalsdta[j]=pvalsdta[j]+sum(TS_data[j]<A[,j])/nrow(A) 
  pvalsdta=c(pvalsdta, outchi)
  names(pvalsdta)=allMethods
  if(TSextra$NoDensity) pvalsdta=pvalsdta[-IndexNotIncluded]
  if(maxProcessor==1) {
    tmp=simpvals(x, TS, typeTS, TSextra, A,  
                 Ranges, nbins, minexpcount, B=B[2])
    pvalsTS=tmp$pvalsTS
    pvalsChi=tmp$pvalsChi
  }
  else {
      z=parallel::clusterCall(cl, simpvals, x, TS, 
                  typeTS, TSextra, A,  Ranges, 
                  nbins, minexpcount, round(B[2]/maxProcessor))
      pvalsTS=z[[1]][[1]]
      pvalsChi=z[[1]][[2]]
      for(i in 2:maxProcessor) {
        pvalsTS=rbind(pvalsTS, z[[i]][[1]])
        pvalsChi=rbind(pvalsChi, z[[i]][[2]])
      }  
      B[1]=nrow(pvalsTS)
      parallel::stopCluster(cl)
  }
  if(TSextra$NoDensity) pvalsTS=pvalsTS[ ,-IndexNotIncluded]
  pvals=cbind(pvalsTS, pvalsChi) 
  pvals=pvals[ ,doMethods,drop=FALSE]
  pvalsdta=pvalsdta[doMethods]
  minp_x=min(pvalsdta)
  minp_sim=apply(pvals[, ,drop=FALSE], 1, min)
  z=seq(0, 1, length=250)
  y=0*z
  for(i in 2:250) y[i]=sum(minp_sim<=z[i])/length(minp_sim)
  I=c(1:250)[z>minp_x][1]-1
  slope=(y[I+1]-y[I])/(z[I+1]-z[I])
  minp_adj=round(y[I]+slope*(minp_x-z[I]),4)
  message("p values of individual tests:")
  for(i in seq_along(pvalsdta)) 
      message(paste(names(pvalsdta)[i],": ", round(pvalsdta[i],4)))
  message(paste0("adjusted p value of combined tests: ",  minp_adj))
  c(pvalsdta, "Min p"=minp_adj)
}
