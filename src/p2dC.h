#ifndef ORDER_H
#define ORDER_H

#include <Rcpp.h>
Rcpp::NumericMatrix order(Rcpp::NumericVector x,
                          Rcpp::NumericVector y);

#endif

#ifndef P2DC_H
#define P2DC_H

#include <Rcpp.h>
Rcpp::NumericMatrix p2dC(Rcpp::NumericMatrix x,
                         Rcpp::Function cdf,
                         Rcpp::NumericVector p,
                         Rcpp::NumericVector Fx);

#endif
