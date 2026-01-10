#' Rearrange 2D discrete data
#' 
#' This function changes a discrete data set given as a nXm counting matrix to a nmX3 matrix
#' 
#' @param x  a matrix of discrete data.
#' @return a rearranged matrix 
#' @export
sq2rec=function(x) {
  xvals=as.numeric(rownames(x))
  yvals=as.numeric(colnames(x))
  cbind(rep(xvals, length(yvals)),
        rep(yvals, each=length(xvals)),
        c(x))
}
