#' Power estimation of tests that find p values
#' 
#' This function finds the power for tests that find their own p values
#' 
#' @param  pnull cdf function
#' @param  ralt function that generates data 
#' @param  param_alt parameters for ralt
#' @param  TS routine that runs the test and returns p values
#' @param  TSextra =list(aa=0), a list of things passed to TS, if needed
#' @param  alpha =0.05 type I error probability of test
#' @param  B =1000 number of simulation runs
#' @keywords internal
#' @return A matrix of power values
#' @export
power_pvals = function(pnull, ralt, param_alt, TS, TSextra=list(aa=0),
                       alpha=0.05, B=1000) {
  
   dta = ralt(param_alt[1])
   typeTS=ifelse(length(formals(TS))==3, 1, 2)
   WithEstimation=ifelse(length(formals(pnull))==1, FALSE, TRUE)
   if(typeTS==1 && (!WithEstimation)) 
      localTS=function(x) TS(x, pnull, 0)
   if(typeTS==2 && (!WithEstimation)) 
     localTS=function(x) TS(x, pnull, 0, TSextra)
   if(typeTS==1 && WithEstimation) 
     localTS=function(x) TS(x, pnull, TSextra$phat(x))
   if(typeTS==2 && WithEstimation) 
     localTS=function(x) TS(x, pnull, TSextra$phat(x), TSextra)
   TS_data=localTS(dta)
   numtests=length(TS_data)

   pwr=matrix(0, length(param_alt), numtests)
   colnames(pwr) = names(TS_data)
   rownames(pwr)=param_alt
   for(i in 1:length(param_alt)) {
     for(j in 1:B) {
         dta=ralt(param_alt[i])
         TS_sim = localTS(dta)
         pwr[i, 1:numtests] = pwr[i, 1:numtests] + ifelse(TS_sim<alpha, 1, 0)
       }
       pwr[i, ] = pwr[i, ,drop=FALSE]/B
  } 
  pwr
}
