#' Power Estimation of Chi Square Tests
#' 
#' This function finds the power of various chi-square tests. 
#' 
#' @param  pnull distribution function to find cdf under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat =function(x) -99, function to estimate parameters 
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Ranges  =matrix(c(-Inf, Inf, -Inf, Inf),2,2), a 2x2 matrix with lower and upper bounds, if any
#' @param  nbins =c(5, 5), number of bins for chi square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  minexpcount =5 minimal expected bin count required
#' @param  dnull =function(x) -99, density function to find probabilities under null hypothesis, mostly used for discrete data, or -99 if missing. 
#' @param  Retry =TRUE, retry if test fails?
#' @param  SuppressMessages =TRUE, should info be shown?
#' @param  B =1000 number of simulation runs to find power
#' @return A numeric matrix of power values.

chi_power = function(
        pnull, ralt, param_alt, phat=function(x) -99, 
        alpha=0.05, Ranges=matrix(c(-Inf, Inf, -Inf, Inf),2,2), 
        nbins=c(5, 5), rate=0, minexpcount=5, dnull=function(x) -99,
        Retry=TRUE, SuppressMessages=TRUE, B=1000) {

   x = ralt(param_alt[1])
   Continuous=TRUE 
   if(ncol(x)>2) {
       if(all(round(x[,3])==x[,3])) Continuous=FALSE
       else {
         message("Chi-square test is only implemented for two dimensional data!")
         return(NULL)
       }
   }
   if(Continuous) {
     pwr=matrix(0, length(param_alt), 2)
     colnames(pwr) = c("ES", "EP")
     rownames(pwr) = param_alt
   }
   else {
     pwr=rep(0, length(param_alt))
     names(pwr) = param_alt
   }
   for(i in 1:length(param_alt)) {
     for(j in 1:B) {
       repeat {
         x = ralt(param_alt[i])
         if(Continuous) {
             tmp = chi_cont_test(x, pnull, phat, Ranges, nbins, 
                                 minexpcount, SuppressMessages)
             if(!is.null(tmp[1])) {
               pwr[i, ] = pwr[i, ] + ifelse(tmp[, 2]<alpha, 1, 0)
               break
             }
             else if(!Retry) return(NULL)
         }     
         else {
             tmp = chi_disc_test(x, pnull, dnull, phat, minexpcount, SuppressMessages)
             if(!is.null(tmp[1])) {
               pwr[i] = pwr[i] + ifelse(tmp[2]<alpha, 1, 0)
               break
             }
             else if(!Retry) return(NULL)
         }     
       }   
     }
   }
   pwr/B
}
