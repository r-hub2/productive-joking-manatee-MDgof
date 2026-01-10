#' Example for a new test 
#' 
#' This shows how a new test routine can be used with Mgof, based on chi square tests
#' 
#' @param x  a data set.
#' @param pnull function to calculate expected counts.
#' @param p parameter for pnull, if needed
#' @param TSextra a list with setup info
#' @return a vector with either values of the test statistic, or p values
#' @export
newTS=function(x, pnull, p, TSextra) {
  repeat {
    if(TSextra$Continuous & (!TSextra$WithEstimation)) 
      tmp=MDgof::chi_cont_test(x, pnull, nbins=c(3, 3),
                               SuppressMessages=TRUE)
    if(TSextra$Continuous & TSextra$WithEstimation) 
      tmp=MDgof::chi_cont_test(x, pnull, phat=TSextra$phat, 
                    nbins=c(3, 3), SuppressMessages=TRUE)
    if((!TSextra$Continuous) & (!TSextra$WithEstimation)) 
      tmp=MDgof::chi_disc_test(x, pnull,
                               SuppressMessages=TRUE)
    if((!TSextra$Continuous) & TSextra$WithEstimation) 
      tmp=MDgof::chi_disc_test(x, pnull, phat=TSextra$phat,
                               SuppressMessages=TRUE)
    if(!is.null(tmp[1])) break
  }
  if(TSextra$Continuous) {
    if(TSextra$Withpvalue) out=tmp[,2]
    else out=tmp[,1]
    names(out)=c("ES33", "EP33")
  }
  else {
    if(TSextra$Withpvalue) out=tmp[2]
    else out=tmp[1]
    names(out)="ES33"
  }
  out
}
