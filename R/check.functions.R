#' Sanity Checks
#' 
#' This function checks whether the inputs have the correct format
#' 
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  phat   =function(x) -99, function to estimate parameters from the data, or -99
#' @param  x      matrix with data
#' @return NULL
#' 
check.functions=function(pnull, rnull, phat=function(x) -99,  x) {
  if(!is.matrix(x)) {
      message("x has to be a matrix")
  }
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  Continuous=TRUE
  if("vals" %in% colnames(x)) {
     Continuous=FALSE
     if(!any(is.wholenumber(x[,colnames(x)!="vals"]))) {
        message("x has a column called vals but some numbers are not counts. Please change the name of that column")
     }
  }  
  npnull = length(formals(pnull))
  nrnull = length(formals(rnull))
# Continuous -  No Estimation 
  if(phat(x)[1]==-99 & Continuous) {
    if(npnull!=1) message("pnull should have one argument for simple hypothesis if data is continuous")
    if(nrnull!=0) message("rnull should have no argument for simple hypothesis")
  } 
# Continuous -  With Estimation 
  if(phat(x)[1]!=-99 & Continuous) {
    if(npnull!=2) message("pnull should have two arguments for composite hypothesis if data is continuous")
    if(nrnull!=1) message("rnull should have one argument for composite hypothesis")
    if(is.function(phat) & length(formals(phat))!=1) 
        message("phat should have one argument x, the data")
  }  
# Discrete -  No Estimation 
  if(phat(x)[1]==-99 && !Continuous) {
    if(npnull!=0) message("pnull should have no argument for simple hypothesis if data is discrete")
    if(nrnull!=0) message("rnull should have no argument for simple hypothesis")
    if(min(diff(pnull()))<0) message("For discrete data pnull should be strictly increasing")
  } 
# Discrete -  With Estimation 
  if(phat(x)[1]!=-99 && !Continuous) {
    if(npnull!=1) message("pnull should have one argument for composite hypothesis if data is discrete")
    if(nrnull!=1) message("rnull should have one argument for composite hypothesis")
    if(is.function(phat) && length(formals(phat))!=1) 
        message("phat should have one argument x, the data")
  }  
  
}
