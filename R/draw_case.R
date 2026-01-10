#' Create plot for any case study
#'
#' This function illustrates any of the case studies.
#' 
#' @param which name or number of the case study.
#' @param Continuous = TRUE for continuous data
#' @param WithEstimation  =FALSE, with parameter estimation
#' @param Dim =2 dimension of data
#' @param palt parameter for alternative. If missing value in study is used.
#' @param nsample =250, sample size.
#' @param Dms =c(1,2, which dimensions are to be shown (for 5D data).
#' @param AltOnly = FALSE show only graph for alternative?
#' @return NULL
#' @export
draw_case=function(which, Continuous=TRUE, 
                   WithEstimation=FALSE, Dim=2, palt, 
                   nsample=1e3, Dms=c(1,2), AltOnly=FALSE) {
  tmp=MDgof::case.studies(which, Continuous, WithEstimation, Dim, nsample)
  x=tmp$rnull()[, Dms]
  if(missing(palt)) {
    palt=tmp$param_alt[2]
    message(paste("Using", palt, "for alternative"))
  }
  y=tmp$ralt(palt)[, Dms]
  z=rbind(x, y)
  dta=data.frame(x=z[,1],
                 y=z[,2],
                 Which=factor(
                   rep(c("Null", "Alt"), each=nsample),
                   ordered = TRUE,
                   levels = c("Null", "Alt")))
  if(AltOnly)
    print(ggplot2::ggplot(dta[(nsample+1):(2*nsample), ], 
                          ggplot2::aes(x=x, y=y))+
            ggplot2::geom_point()+ggplot2::ggtitle("Alt"))
  else  print(ggplot2::ggplot(dta, ggplot2::aes(x=x, y=y))+
                ggplot2::geom_point()+ggplot2::facet_grid(.~Which))
  
} 
