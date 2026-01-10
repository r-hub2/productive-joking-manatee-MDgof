#' Create various case studies
#'
#' This function creates the functions needed to run the various case studies.
#' 
#' @param which name or number of the case study.
#' @param Continuous = TRUE for continuous data
#' @param WithEstimation  =FALSE, with parameter estimation
#' @param Dim =2 dimension of data
#' @param nsample =250, sample size.
#' @param nbins = c(5,5), number of bins for discretized data
#' @param bintype = "ES", ES for equal spaced bins or EP for equal probability bins.
#' @param ReturnCaseNames =FALSE, just return names of case studies?
#' @return a list of functions
#' @export
case.studies=function(which, Continuous=TRUE, WithEstimation=FALSE, Dim=2,
                      nsample=250, nbins=c(5,5), bintype = "ES",
                      ReturnCaseNames=FALSE) {

  if(ReturnCaseNames) {
    if(Dim==5)
      return(MDgof::case.studies.cont.D5(ReturnCaseNames=TRUE))    
    if(!WithEstimation)
      return(MDgof::case.studies.cont(ReturnCaseNames=TRUE))
    return(MDgof::case.studies.est(ReturnCaseNames=TRUE))    
   }
   if(Dim==5) {
     return(MDgof::case.studies.cont.D5(which, nsample=nsample))
   }
   if(Continuous & (!WithEstimation)) {
      return(MDgof::case.studies.cont(which, nsample=nsample))
   }
   if(Continuous & WithEstimation) {
     return(MDgof::case.studies.est(which, nsample=nsample))
   }
   if(!Continuous) {
     return(MDgof::case.studies.disc(which, WithEstimation, bintype,
                              nbins, nsample=nsample))
   }
  
}

#' Create copula objects
#'
#' This function creates copula objects
#' 
#' @param family name of copula.
#' @param p parameter of copula.
#' @param d dimension
#' @return a copula object
#' @export
gen.cop=function(family, p, d=2) {
  if(family=="uniform") 
    return(copula::indepCopula(dim=d))
  if(family=="HuslerReiss" | family=="MarshallOlkin") {
    if(family=="HuslerReiss")
      copula_fun = utils::getFromNamespace("huslerReissCopula", "copula")
    if(family=="MarshallOlkin")  
      copula_fun = utils::getFromNamespace("moCopula", "copula")
  }
  else    
    copula_fun <- utils::getFromNamespace(paste0(tolower(family), "Copula"), "copula")
  if(missing(d)) return(copula_fun(p))
  return(copula_fun(p, dim=d))
}

