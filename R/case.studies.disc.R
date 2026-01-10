# create a new environment to store the parameter estimation from the continuous case study.
phat.env=new.env(parent = emptyenv())

#' Discretize 2D data from case studies
#' 
#' This function provides the info necessary to run the case studies for discrete data.
#' 
#' @param which name or number of desired case study.
#' @param WithEstimation = FALSE, case study with or without parameter estimation.
#' @param bintype = "ES", ES for equal spaced bins or EP for equal probability bins.
#' @param nbins =c(5, 5) number of bins to use in x and y direction 
#' @param nsample = 250, required sample size
#' @return a list with needed stuff
#' @export
case.studies.disc=function(which, WithEstimation=FALSE, 
                           bintype = "ES",  nbins=c(5, 5), nsample=250) {
  
  if(WithEstimation) 
     list.of.cases=MDgof::case.studies.est(ReturnCaseNames=TRUE)
  else
    list.of.cases=MDgof::case.studies(ReturnCaseNames=TRUE)
  if(is.numeric(which[1])) which=list.of.cases[which]
  if(!WithEstimation) {
      tmp=MDgof::case.studies(which, nsample=nsample) 
      rnull1=function() 
         MDgof::discretize(tmp$rnull(), tmp$Range, nbins=nbins)
      ralt=function(p) 
        MDgof::discretize(tmp$ralt(p), tmp$Range,  nbins=nbins)
      return(list(pnull=tmp$pnull, rnull=rnull1, ralt=ralt,
                  param_alt=tmp$param_alt))
  }
  else {
     tmp=MDgof::case.studies.est(which, nsample=nsample)
     rnull2=function(p) {
       x=tmp$rnull(p)
       phat.env$phat=tmp$phat(x)
       MDgof::discretize(x, tmp$Range, nbins=nbins)
     } 
     ralt=function(p) {
       x=tmp$ralt(p)
       phat.env$phat=tmp$phat(x)
       MDgof::discretize(x, tmp$Range, nbins=nbins)
     } 
     phat=function(x) return(phat.env$phat)
     return(list(pnull=tmp$pnull, rnull=rnull2, ralt=ralt,
                 phat=phat, param_alt=tmp$param_alt))
  }
  
}
