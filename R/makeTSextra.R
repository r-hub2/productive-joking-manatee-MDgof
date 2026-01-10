#' Create list with needed info
#' 
#' This function creates a list with info  needed in various parts of the package
#' 
#' @param  x data set
#' @param  Continuous =TRUE, is data continuous?
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  phat  =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  dnull =function(x) -99, density function under the null hypothesis, if available, or -99 if missing
#' @param  Ranges Range of variables
#' @param  TSextra (optional) list passed to TS, if needed.
#' @return A list with vectors of test statistics and p.values
#' @export
makeTSextra <- function(x, Continuous, pnull, rnull, phat=function(x) -99, 
                    dnull=function(x) -99, Ranges, TSextra) {
  if(missing(TSextra)) {
    TSextra=list(pnull=pnull, rnull=rnull, phat=phat, Range=Ranges,
                 dnull=dnull, ripleyK=RipleyK, Continuous=Continuous)
    if(Continuous) {
      Eval=gen_eval(rnull, phat(x), 100) #Needed for Bakshaev_Rudzkis
      TSextra=c(TSextra, list(Eval=Eval))
    }   
  }   
  else {
    if(!is.list(TSextra)) {
      message("If TSextra is provided it has to be a list object!")
      return(NULL)
    } 
    TSextra = c(TSextra, pnull=pnull, rnull=rnull, phat=phat, 
                dnull=dnull, ripleyK=RipleyK, 
                Range=Ranges, Continuous=Continuous)
    if(Continuous) {
      Eval=gen_eval(rnull, phat(x), 100) #Needed for Bakshaev_Rudzkis
      TSextra=c(TSextra, list(Eval=Eval))
    }
  }  
  NoDensity=FALSE
  if(length(formals(dnull))==1 && dnull(x)[1]<0) NoDensity=TRUE
  if(ncol(x)>2) NoDensity=TRUE
  TSextra=c(TSextra, NoDensity=NoDensity)
  TSextra$knn=function(x) -99
  if(!NoDensity) {
    if(Continuous) {
      TSextra$knn=function(x) FNN::get.knn(x, 1)
      TSextra$nn=FNN::get.knn(x, 1)
      TSextra$scf=function(x, p=0) {
        if(length(formals(dnull))==1) f=function(x) dnull(x)
        if(length(formals(dnull))==2) f=function(x) dnull(x, p)
        if(!is.matrix(x)) return(grad_vec(x, f))
        grad_mat(x, f)
      }
    }
  }
  return(TSextra)
  
}
