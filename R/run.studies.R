#' Benchmarking for Multivariate Goodness-of-fit Tests
#' 
#' This function runs the case studies included in the package.
#' 
#' For details consult vignette(package="MDgof")
#' 
#' @param study either the name of the study, or its number in the list. If missing all the studies are run.
#' @param Continuous =TRUE, run cases for continuous data.
#' @param WithEstimation =FALSE, run case studies with or without parameter estimation?
#' @param Dim =2 two or five-dimensional continuous data sets?
#' @param TS routine to calculate new test statistics. 
#' @param TSextra list passed to TS (optional).
#' @param With.p.value =FALSE, does user supplied routine return p values?
#' @param nsample = 250, desired sample size. 250 is used in included case studies.
#' @param nbins =c(5,5) number of bins for discretized data.
#' @param alpha =0.05,  type I error probability of tests. 0.05 is used in included case studies.
#' @param param_alt (list of) values of parameter under the alternative hypothesis. 
#'                  If missing included values are used.
#' @param SuppressMessages =TRUE, should informative messages be printed?
#' @param B = 1000, number of simulation runs.
#' @param  maxProcessor  number of cores to use. If missing the number of physical cores-1 
#'             is used. If set to 1 no parallel processing is done.
#' @return A (list of ) matrices of p.values.
#' @examples
#' #Examples are run with a super small B=25 simulation runs to satisfy CRAN submission rules.
#' #Run a new test for studies 1-3 for continuous data and without estimation.
#' #The new test is an (included) chi square test that finds it's own p value.
#' TSextra=list(Continuous=TRUE, WithEstimation=FALSE, Withpvalue=TRUE)
#' MDgof::run.studies(Continuous=TRUE, WithEstimation=FALSE, 
#'            study=1:3, TS=MDgof::newTS, TSextra=TSextra, 
#'            With.p.value = TRUE, B=25, maxProcessor = 1)
#' #Run included tests for studies 1-3 for discrete data and without estimation,
#' #but with type I error alpha=0.1
#' p=MDgof::power_studies_disc_results[[3]][1:3,,drop=FALSE]    
#' MDgof::run.studies(Continuous=FALSE, WithEstimation=FALSE, 
#'            study=1:3, param_alt=p,alpha=0.1, B=25, maxProcessor = 1)         
#' @export
run.studies <- function(study, Continuous=TRUE, WithEstimation=FALSE, 
                      Dim=2, TS, TSextra, With.p.value=FALSE,  
                      nsample=250, nbins=c(5,5), 
                      alpha=0.05, param_alt, 
                      SuppressMessages =TRUE, B=1000, maxProcessor) {
  list.of.studies=MDgof::case.studies(WithEstimation=WithEstimation,
                                      Dim=Dim, ReturnCaseNames=TRUE)
  NewTS=ifelse(missing(TS), FALSE, TRUE)
  if(missing(study)) study=1:length(list.of.studies) #Do all of them
  if(is.numeric(study)) study=list.of.studies[study] #get study names
  require(MDgof)
  pwr=get(paste0("power_studies",
                 ifelse(Continuous,"_cont", "_disc"),
                 ifelse(WithEstimation, "_est", ""),
                 ifelse(Dim==5, "_D5", ""),
                 "_results"), envir = as.environment("package:MDgof"))
  if(missing(param_alt)) {
     param_alt=pwr[[3]][study, ,drop=FALSE]
  }
  else {
     if(length(study)>1) {
        if(is.matrix(param_alt)) {
           if(nrow(param_alt)!=length(study)) {
              message("param_alt should be a matrix with one row for each study")
              return(NULL)
           }   
        }
        else {
          message("if more than one study is run param_alt has to be a matrix with the one row for each study")
          return(NULL)
        }  
      }
      else param_alt=rbind(param_alt)
  } 
  if(NewTS) {
    param_alt=param_alt[,2]
    allpwr=pwr[[1]][study, , drop=FALSE]
  }
  phat=function(x) -99
  for(i in seq_along(study)) {
        tmp=MDgof::case.studies(study[i], Continuous, WithEstimation, Dim,
                          nsample=nsample, nbins=nbins)
        if(WithEstimation) {
          phat=tmp$phat
          if(!NewTS) param_alt[i,1]=-1
        }
        if(Continuous) {
          Ra=tmp$Range
          dnull=tmp$dnull
        }  
        else {Ra=NA; dnull=function(x) -99}
        if(!SuppressMessages) {
            message(paste("Running case study", study[i], 
                          "for", ifelse(Continuous, "continuous", "discrete"),
                          " ", ifelse(Dim==2, 2, 5),"dimensional",
                          "data and ", ifelse(WithEstimation, "with", "without"),
                          "parameter estimation"))
        }
        if(!NewTS) pr=param_alt[i, ]
        if(NewTS) pr=param_alt[i]
        newpwr=MDgof::gof_power(tmp$pnull, tmp$rnull, tmp$ralt, param_alt=pr, 
                             alpha=alpha,  dnull=dnull, TS=TS, TSextra=TSextra,
                             Ranges=Ra, phat=phat, With.p.value=With.p.value,
                             SuppressMessages=SuppressMessages, 
                             nbins=nbins, B=B, maxProcessor=maxProcessor) 
        if(!NewTS) {
          pwr[["Null"]][study[i], ]=newpwr[1, ]
          pwr[["Pow"]][study[i], ]=newpwr[2, ]
        }
        else {
          if(i==1)  {
            allpwr=cbind(matrix(0, length(study), ncol(newpwr)), allpwr)
            colnames(allpwr)[1:ncol(newpwr)] = colnames(newpwr)
          }  
          allpwr[study[i], 1:ncol(newpwr)]=newpwr
        }
  }  
  if(!NewTS) return(pwr)
  if(length(study)>1) {
    a1=apply(allpwr, 1, rank)
    message("Average number of times a test is close to best:")
    print(sort(apply(a1,1,mean)))
  }  
  allpwr
}
