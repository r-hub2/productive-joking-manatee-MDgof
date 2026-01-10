#ifndef ROSENBLATTC_H
#define ROSENBLATTC_H

#include <Rcpp.h>
Rcpp::NumericMatrix rosenblattC(
       Rcpp::NumericMatrix x, 
       Rcpp::Function cdf, 
       Rcpp::NumericVector p,
       Rcpp::NumericMatrix Range
        );

#endif
